#!/usr/bin/env python3
"""
基盤地図情報 DEM XML → STL 変換スクリプト
==========================================
国土地理院の基盤地図情報DEMデータ（XML形式）を
3Dプリント可能なSTLファイルに変換します。

対応データ:
  - DEM5A（航空レーザ測量, 5mメッシュ, 150×225グリッド）
  - DEM5B/5C（写真測量, 5mメッシュ）
  - DEM10B（10mメッシュ, 750×1125グリッド）
  - 上記の混在データ（自動リサンプルで解像度を統一）

機能:
  - GeoJSON行政境界による切り抜き（任意の市区町村）
  - 市域形状のみのSTL出力（余白なし）または四角い台座付き
  - 5m/10mメッシュの混在を自動リサンプルで統一
  - ZIPファイル（入れ子含む）の自動展開
  - メモリ制限付き自動ダウンサンプル

依存: numpy のみ（pip不要の標準的な環境で動作）

使い方:
  python fgd_dem_to_stl.py <DEMフォルダ or ZIP> \\
    -g <行政境界GeoJSON> \\
    -n <市区町村名> \\
    -o output.stl \\
    -w 160 -z 2.0 --max-points 10000000

  # 四角い台座付きバージョン
  python fgd_dem_to_stl.py <DEMフォルダ> -g boundary.geojson -n 糸魚川市 --rect

  # 境界なし（全データを四角いSTLに）
  python fgd_dem_to_stl.py <DEMフォルダ> -o terrain.stl
"""

import argparse
import glob
import json
import os
import struct
import sys
import zipfile
import tempfile
from xml.etree import ElementTree as ET
import numpy as np


# =====================================================================
#  GeoJSON 境界処理
# =====================================================================

def load_boundary(geojson_path, area_name):
    """
    GeoJSONから指定地域のポリゴンを抽出する。

    Args:
        geojson_path: GeoJSONファイルパス
        area_name: 検索する地域名（例: "糸魚川", "富士見町"）
                   propertiesの値のいずれかに部分一致すればヒット
    Returns:
        polygons: [(outer_ring, [hole_rings]), ...] のリスト
        bbox: {'lat_min', 'lat_max', 'lon_min', 'lon_max'}
    """
    with open(geojson_path, 'r', encoding='utf-8') as f:
        data = json.load(f)

    polygons = []
    for feat in data['features']:
        props = feat.get('properties', {})
        is_target = any(v and area_name in str(v) for v in props.values())
        if not is_target:
            continue

        geom = feat['geometry']
        if geom['type'] == 'Polygon':
            outer = np.array(geom['coordinates'][0], dtype=np.float64)
            holes = [np.array(h, dtype=np.float64) for h in geom['coordinates'][1:]]
            polygons.append((outer, holes))
        elif geom['type'] == 'MultiPolygon':
            for poly_coords in geom['coordinates']:
                outer = np.array(poly_coords[0], dtype=np.float64)
                holes = [np.array(h, dtype=np.float64) for h in poly_coords[1:]]
                polygons.append((outer, holes))

    if not polygons:
        print(f"エラー: GeoJSON内に「{area_name}」を含むフィーチャーが見つかりません。")
        print("  GeoJSON内の地域名を確認してください。")
        # デバッグ: 最初の5件のproperties表示
        for i, feat in enumerate(data['features'][:5]):
            print(f"  例 [{i}]: {feat.get('properties', {})}")
        sys.exit(1)

    print(f"「{area_name}」ポリゴン: {len(polygons)} 個")

    all_lons = np.concatenate([p[0][:, 0] for p in polygons])
    all_lats = np.concatenate([p[0][:, 1] for p in polygons])
    bbox = {
        'lat_min': all_lats.min(), 'lat_max': all_lats.max(),
        'lon_min': all_lons.min(), 'lon_max': all_lons.max(),
    }
    print(f"境界範囲: 緯度 {bbox['lat_min']:.4f}～{bbox['lat_max']:.4f}, "
          f"経度 {bbox['lon_min']:.4f}～{bbox['lon_max']:.4f}")

    return polygons, bbox


def points_in_polygon(points_lon, points_lat, polygon):
    """Ray casting法（numpy ベクトル化版）"""
    n = len(polygon)
    inside = np.zeros(len(points_lon), dtype=bool)
    j = n - 1
    for i in range(n):
        xi, yi = polygon[i, 0], polygon[i, 1]
        xj, yj = polygon[j, 0], polygon[j, 1]
        cond1 = (yi > points_lat) != (yj > points_lat)
        if cond1.any():
            slope = (xj - xi) / (yj - yi + 1e-30)
            x_intersect = xi + slope * (points_lat - yi)
            cross = cond1 & (points_lon < x_intersect)
            inside ^= cross
        j = i
    return inside


def create_boundary_mask(grid_lats, grid_lons, polygons):
    """グリッド上の各点が境界ポリゴン内にあるかのマスクを生成"""
    nrows = len(grid_lats)
    ncols = len(grid_lons)
    lon_grid, lat_grid = np.meshgrid(grid_lons, grid_lats)
    flat_lon = lon_grid.flatten()
    flat_lat = lat_grid.flatten()

    print(f"マスク計算中... ({nrows}x{ncols} = {len(flat_lon):,} points)")
    mask = np.zeros(len(flat_lon), dtype=bool)

    for i, (outer, holes) in enumerate(polygons):
        if i % 10 == 0:
            print(f"  ポリゴン {i}/{len(polygons)}...")
        olon_min, olon_max = outer[:, 0].min(), outer[:, 0].max()
        olat_min, olat_max = outer[:, 1].min(), outer[:, 1].max()
        bbox_mask = ((flat_lon >= olon_min) & (flat_lon <= olon_max) &
                     (flat_lat >= olat_min) & (flat_lat <= olat_max))
        if not bbox_mask.any():
            continue
        subset_lon = flat_lon[bbox_mask]
        subset_lat = flat_lat[bbox_mask]
        in_outer = points_in_polygon(subset_lon, subset_lat, outer)
        for hole in holes:
            in_hole = points_in_polygon(subset_lon, subset_lat, hole)
            in_outer &= ~in_hole
        bbox_indices = np.where(bbox_mask)[0]
        mask[bbox_indices] |= in_outer

    mask_2d = mask.reshape(nrows, ncols)
    print(f"境界内ポイント: {mask_2d.sum():,}/{nrows*ncols:,} "
          f"({mask_2d.sum()/(nrows*ncols)*100:.1f}%)")
    return mask_2d


# =====================================================================
#  DEM XML パーサー
# =====================================================================

def parse_fgd_dem_xml(filepath):
    """
    基盤地図情報DEM XMLファイルをパースして標高データを返す。
    DEM5A/5B/5C (5mメッシュ) と DEM10B (10mメッシュ) に対応。
    """
    try:
        tree = ET.parse(filepath)
    except ET.ParseError:
        return None

    root = tree.getroot()
    ns = {'gml': 'http://www.opengis.net/gml/3.2',
          'fgd': 'http://fgd.gsi.go.jp/spec/2008/FGD_GMLSchema'}

    def find_text(paths):
        for p in paths:
            e = root.find(p, ns)
            if e is not None and e.text:
                return e.text.strip()
        return None

    mesh = find_text(['.//fgd:mesh',
                      './/{http://fgd.gsi.go.jp/spec/2008/FGD_GMLSchema}mesh'])
    if not mesh:
        return None

    lower = find_text(['.//gml:lowerCorner'])
    upper = find_text(['.//gml:upperCorner'])
    high = find_text(['.//gml:high'])
    tuples = find_text(['.//gml:tupleList'])

    if not all([lower, upper, high, tuples]):
        return None

    lp = lower.split()
    up = upper.split()
    hp = high.split()

    lat_min, lon_min = float(lp[0]), float(lp[1])
    lat_max, lon_max = float(up[0]), float(up[1])
    ncols = int(hp[0]) + 1
    nrows = int(hp[1]) + 1

    elevations = np.full(nrows * ncols, np.nan, dtype=np.float32)
    idx = 0
    for line in tuples.split('\n'):
        line = line.strip()
        if not line:
            continue
        comma_pos = line.rfind(',')
        if comma_pos < 0:
            idx += 1
            continue
        try:
            val = float(line[comma_pos+1:])
            if val > -9000:
                elevations[idx] = val
        except ValueError:
            pass
        idx += 1
        if idx >= nrows * ncols:
            break

    return {
        'mesh': mesh,
        'lat_min': lat_min, 'lat_max': lat_max,
        'lon_min': lon_min, 'lon_max': lon_max,
        'nrows': nrows, 'ncols': ncols,
        'elevations': elevations
    }


# =====================================================================
#  バイリニア リサンプル
# =====================================================================

def resample_tile_fast(tile_2d, src_nrows, src_ncols, dst_nrows, dst_ncols):
    """
    numpy ベクトル化版の bilinear リサンプル。
    10mメッシュタイルを5mベースのグリッドに統一するために使用。
    NaN（欠損値）を考慮した重み付き補間。
    """
    if src_nrows == dst_nrows and src_ncols == dst_ncols:
        return tile_2d

    src_r = np.linspace(0, src_nrows - 1, dst_nrows).astype(np.float32)
    src_c = np.linspace(0, src_ncols - 1, dst_ncols).astype(np.float32)

    r0 = np.floor(src_r).astype(int)
    r1 = np.minimum(r0 + 1, src_nrows - 1)
    fr = (src_r - r0).reshape(-1, 1)

    c0 = np.floor(src_c).astype(int)
    c1 = np.minimum(c0 + 1, src_ncols - 1)
    fc = (src_c - c0).reshape(1, -1)

    v00 = tile_2d[np.ix_(r0, c0)]
    v01 = tile_2d[np.ix_(r0, c1)]
    v10 = tile_2d[np.ix_(r1, c0)]
    v11 = tile_2d[np.ix_(r1, c1)]

    w00 = (1 - fr) * (1 - fc)
    w01 = (1 - fr) * fc
    w10 = fr * (1 - fc)
    w11 = fr * fc

    valid00 = ~np.isnan(v00)
    valid01 = ~np.isnan(v01)
    valid10 = ~np.isnan(v10)
    valid11 = ~np.isnan(v11)

    sum_w = (w00 * valid00 + w01 * valid01 + w10 * valid10 + w11 * valid11)
    sum_v = (w00 * np.where(valid00, v00, 0) +
             w01 * np.where(valid01, v01, 0) +
             w10 * np.where(valid10, v10, 0) +
             w11 * np.where(valid11, v11, 0))

    with np.errstate(invalid='ignore'):
        result = np.where(sum_w > 0, sum_v / sum_w, np.nan)

    return result.astype(np.float32)


# =====================================================================
#  ファイル展開
# =====================================================================

def extract_all_xmls(input_path):
    """DEMフォルダまたはZIPからXMLファイル一覧を取得（入れ子ZIP対応）"""
    if os.path.isdir(input_path):
        xmls = glob.glob(os.path.join(input_path, '**', '*.xml'), recursive=True)
        return [x for x in xmls
                if 'DEM' in os.path.basename(x).upper()
                or 'dem' in os.path.basename(x)]
    if zipfile.is_zipfile(input_path):
        tmpdir = tempfile.mkdtemp()
        with zipfile.ZipFile(input_path, 'r') as zf:
            zf.extractall(tmpdir)
        # 入れ子ZIPを展開
        for root_d, dirs, files in os.walk(tmpdir):
            for f in files:
                fp = os.path.join(root_d, f)
                if f.endswith('.zip'):
                    try:
                        with zipfile.ZipFile(fp, 'r') as zf2:
                            zf2.extractall(tmpdir)
                    except Exception:
                        pass
        xmls = glob.glob(os.path.join(tmpdir, '**', '*.xml'), recursive=True)
        return [x for x in xmls if 'DEM' in os.path.basename(x).upper()]
    return [input_path] if input_path.endswith('.xml') else []


# =====================================================================
#  グリッド構築
# =====================================================================

def build_grid(xml_files, bbox=None, max_points=4_000_000):
    """
    DEMタイルを読み込んで統合グリッドを構築する。

    Args:
        xml_files: XMLファイルパスのリスト
        bbox: 切り抜き範囲 {'lat_min','lat_max','lon_min','lon_max'} or None（全体）
        max_points: グリッドの最大ポイント数（メモリ制限）
    """
    print(f"\nXMLファイル数: {len(xml_files)}")

    # 切り抜き範囲
    if bbox:
        margin = 0.002
        clip_lat_min = bbox['lat_min'] - margin
        clip_lat_max = bbox['lat_max'] + margin
        clip_lon_min = bbox['lon_min'] - margin
        clip_lon_max = bbox['lon_max'] + margin
    else:
        clip_lat_min = clip_lat_max = clip_lon_min = clip_lon_max = None

    # Pass 1: メタデータ収集
    print("Pass 1: タイルを選別中...")
    tile_meta = []
    skipped = 0
    for i, f in enumerate(sorted(xml_files)):
        if i % 200 == 0:
            print(f"  {i}/{len(xml_files)}...")
        t = parse_fgd_dem_xml(f)
        if t is None:
            continue
        if bbox:
            if (t['lat_max'] < clip_lat_min or t['lat_min'] > clip_lat_max or
                t['lon_max'] < clip_lon_min or t['lon_min'] > clip_lon_max):
                skipped += 1
                continue
        tile_meta.append((f, t['lat_min'], t['lat_max'], t['lon_min'], t['lon_max'],
                         t['nrows'], t['ncols']))

    print(f"  対象タイル: {len(tile_meta)} (スキップ: {skipped})")
    if not tile_meta:
        print("エラー: 対象範囲内にタイルがありません")
        sys.exit(1)

    # bboxがない場合、全タイルの範囲を使用
    if not bbox:
        clip_lat_min = min(t[1] for t in tile_meta) - 0.001
        clip_lat_max = max(t[2] for t in tile_meta) + 0.001
        clip_lon_min = min(t[3] for t in tile_meta) - 0.001
        clip_lon_max = max(t[4] for t in tile_meta) + 0.001

    # 出力グリッド解像度（5mタイル基準）
    ref_tile = None
    for tm in tile_meta:
        if tm[5] == 150 and tm[6] == 225:  # DEM5A/5B/5C
            ref_tile = tm
            break
    if ref_tile is None:
        ref_tile = tile_meta[0]

    res_lat = (ref_tile[2] - ref_tile[1]) / ref_tile[5]
    res_lon = (ref_tile[4] - ref_tile[3]) / ref_tile[6]

    total_nrows = int(round((clip_lat_max - clip_lat_min) / res_lat)) + 1
    total_ncols = int(round((clip_lon_max - clip_lon_min) / res_lon)) + 1

    # 自動ダウンサンプル
    downsample = 1
    while (total_nrows // downsample) * (total_ncols // downsample) > max_points:
        downsample += 1

    grid_nrows = total_nrows // downsample + 1
    grid_ncols = total_ncols // downsample + 1
    effective_res_lat = (clip_lat_max - clip_lat_min) / (grid_nrows - 1)
    effective_res_lon = (clip_lon_max - clip_lon_min) / (grid_ncols - 1)

    print(f"\nグリッドサイズ:")
    print(f"  フル解像度: {total_nrows} x {total_ncols}")
    if downsample > 1:
        print(f"  ダウンサンプル: {downsample}倍 → {grid_nrows} x {grid_ncols}")
    print(f"  出力解像度: lat {effective_res_lat:.7f}°, lon {effective_res_lon:.7f}°")

    grid = np.full((grid_nrows, grid_ncols), np.nan, dtype=np.float32)

    # Pass 2: タイル配置（地理座標ベースのリサンプル）
    print("\nPass 2: タイルを配置中...")
    for i, (fpath, lat_min, lat_max, lon_min, lon_max, nrows, ncols) in enumerate(tile_meta):
        if i % 50 == 0:
            print(f"  {i}/{len(tile_meta)}...")
        t = parse_fgd_dem_xml(fpath)
        if t is None:
            continue

        elev = t['elevations']
        tile_2d = elev.reshape(nrows, ncols) if len(elev) == nrows * ncols else None
        if tile_2d is None:
            tile_2d = np.full((nrows, ncols), np.nan, dtype=np.float32)
            tile_2d.flat[:len(elev)] = elev

        # 地理座標ベースでの目標セル数を計算
        dst_nrows = max(1, int(round((lat_max - lat_min) / effective_res_lat)))
        dst_ncols = max(1, int(round((lon_max - lon_min) / effective_res_lon)))

        # サイズが違う場合リサンプル（10m→5m変換含む）
        if nrows != dst_nrows or ncols != dst_ncols:
            tile_2d = resample_tile_fast(tile_2d, nrows, ncols, dst_nrows, dst_ncols)

        # 出力グリッドでの配置位置
        row_off = int(round((clip_lat_max - lat_max) / effective_res_lat))
        col_off = int(round((lon_min - clip_lon_min) / effective_res_lon))

        tr, tc = tile_2d.shape
        r_end = min(row_off + tr, grid_nrows)
        c_end = min(col_off + tc, grid_ncols)
        r_start = max(row_off, 0)
        c_start = max(col_off, 0)

        tile_r_start = r_start - row_off
        tile_c_start = c_start - col_off
        tile_r_end = tile_r_start + (r_end - r_start)
        tile_c_end = tile_c_start + (c_end - c_start)

        if r_end > r_start and c_end > c_start:
            patch = tile_2d[tile_r_start:tile_r_end, tile_c_start:tile_c_end]
            valid = ~np.isnan(patch)
            grid[r_start:r_end, c_start:c_end] = np.where(
                valid, patch, grid[r_start:r_end, c_start:c_end]
            )

    # NaN補間（8方向反復）
    nan_count = np.isnan(grid).sum()
    valid_count = grid.size - nan_count
    print(f"\n有効データ: {valid_count:,}/{grid.size:,} ({valid_count/grid.size*100:.1f}%)")

    if nan_count > 0 and valid_count > 0:
        print("欠損値を補間中（8方向）...")
        for it in range(50):
            nans = np.isnan(grid)
            if not nans.any():
                break
            shifted = np.stack([
                np.roll(grid, 1, axis=0),
                np.roll(grid, -1, axis=0),
                np.roll(grid, 1, axis=1),
                np.roll(grid, -1, axis=1),
                np.roll(np.roll(grid, 1, axis=0), 1, axis=1),
                np.roll(np.roll(grid, 1, axis=0), -1, axis=1),
                np.roll(np.roll(grid, -1, axis=0), 1, axis=1),
                np.roll(np.roll(grid, -1, axis=0), -1, axis=1),
            ])
            with np.errstate(invalid='ignore'):
                counts = np.sum(~np.isnan(shifted), axis=0)
                sums = np.nansum(shifted, axis=0)
                avg = np.where(counts > 0, sums / counts, np.nan)
            grid = np.where(nans, avg, grid)
            remaining = np.isnan(grid).sum()
            if it % 10 == 0:
                print(f"  反復{it+1}: 残り {remaining:,}")
            if remaining == 0:
                break
        if np.isnan(grid).any():
            grid = np.where(np.isnan(grid), np.nanmin(grid), grid)

    # ガウシアンブラー（タイル継ぎ目のなじませ）
    print("\n全体スムージング（ガウシアン 3x3）...")
    kernel = np.array([[1, 2, 1],
                       [2, 4, 2],
                       [1, 2, 1]], dtype=np.float32) / 16.0
    nrows_g, ncols_g = grid.shape
    padded = np.pad(grid, 1, mode='edge')
    result = np.zeros_like(grid)
    for dr in range(3):
        for dc in range(3):
            result += kernel[dr, dc] * padded[dr:dr+nrows_g, dc:dc+ncols_g]
    grid = result

    grid_lats = np.linspace(clip_lat_max, clip_lat_min, grid_nrows)
    grid_lons = np.linspace(clip_lon_min, clip_lon_max, grid_ncols)

    return grid, clip_lat_min, clip_lat_max, clip_lon_min, clip_lon_max, grid_lats, grid_lons


# =====================================================================
#  STL 書き出し
# =====================================================================

def _stl_helpers(f):
    """STLバイナリ書き出し用のヘルパー関数を返す"""
    attr = struct.pack('<H', 0)
    counter = [0]

    def write_tri(n, v1, v2, v3):
        f.write(struct.pack('<3f', *n))
        f.write(struct.pack('<3f', *v1))
        f.write(struct.pack('<3f', *v2))
        f.write(struct.pack('<3f', *v3))
        f.write(attr)
        counter[0] += 1

    def calc_normal(v1, v2, v3):
        ux, uy, uz = v2[0]-v1[0], v2[1]-v1[1], v2[2]-v1[2]
        vx, vy, vz = v3[0]-v1[0], v3[1]-v1[1], v3[2]-v1[2]
        nx = uy*vz - uz*vy
        ny = uz*vx - ux*vz
        nz = ux*vy - uy*vx
        nl = (nx*nx + ny*ny + nz*nz) ** 0.5
        if nl > 0:
            return (nx/nl, ny/nl, nz/nl)
        return (0.0, 0.0, 1.0)

    return write_tri, calc_normal, counter


def grid_to_stl_shape(grid, lat_min, lat_max, lon_min, lon_max,
                      boundary_mask, output_path,
                      model_width_mm=160.0, z_exaggeration=2.0,
                      base_height_mm=3.0):
    """境界形状のみのSTLを生成（余白なし版）"""
    nrows, ncols = grid.shape
    print(f"\nSTL生成中（形状のみ）... ({nrows} x {ncols})")

    mid_lat = (lat_min + lat_max) / 2.0
    lat_span_m = (lat_max - lat_min) * 111320.0
    lon_span_m = (lon_max - lon_min) * 111320.0 * np.cos(np.radians(mid_lat))

    scale = model_width_mm / lon_span_m
    model_height_mm = lat_span_m * scale

    min_elev = float(np.min(grid[boundary_mask])) if boundary_mask.any() else 0.0
    max_elev = float(np.max(grid[boundary_mask])) if boundary_mask.any() else 100.0
    elev_range_mm = (max_elev - min_elev) * scale * z_exaggeration

    print(f"実寸: {lon_span_m/1000:.1f}km x {lat_span_m/1000:.1f}km")
    print(f"モデル: {model_width_mm:.1f} x {model_height_mm:.1f} x "
          f"{elev_range_mm + base_height_mm:.1f} mm")
    print(f"標高: {min_elev:.0f}m ～ {max_elev:.0f}m")

    x = np.linspace(0, model_width_mm, ncols, dtype=np.float32)
    y = np.linspace(model_height_mm, 0, nrows, dtype=np.float32)
    xx, yy = np.meshgrid(x, y)
    zz = ((grid - min_elev) * scale * z_exaggeration + base_height_mm).astype(np.float32)
    zz[~boundary_mask] = base_height_mm

    z_bottom = np.float32(0.0)

    cell_mask = (boundary_mask[:-1, :-1] | boundary_mask[:-1, 1:] |
                 boundary_mask[1:, :-1] | boundary_mask[1:, 1:])
    active_cells = int(cell_mask.sum())
    print(f"アクティブセル: {active_cells:,}/{(nrows-1)*(ncols-1):,}")

    with open(output_path, 'wb') as f:
        header = b'FGD DEM to STL - Shape Mode' + b'\0' * (80 - 27)
        f.write(header)
        count_pos = f.tell()
        f.write(struct.pack('<I', 0))

        write_tri, calc_normal, counter = _stl_helpers(f)

        # 上面
        for r in range(nrows - 1):
            if r % 300 == 0:
                print(f"  上面: {r}/{nrows-1}")
            for c in range(ncols - 1):
                if not cell_mask[r, c]:
                    continue
                v00 = (xx[r,c], yy[r,c], zz[r,c])
                v10 = (xx[r,c+1], yy[r,c+1], zz[r,c+1])
                v01 = (xx[r+1,c], yy[r+1,c], zz[r+1,c])
                v11 = (xx[r+1,c+1], yy[r+1,c+1], zz[r+1,c+1])
                n1 = calc_normal(v00, v10, v01)
                write_tri(n1, v00, v10, v01)
                n2 = calc_normal(v10, v11, v01)
                write_tri(n2, v10, v11, v01)

        # 底面
        print("  底面...")
        for r in range(nrows - 1):
            for c in range(ncols - 1):
                if not cell_mask[r, c]:
                    continue
                v00 = (xx[r,c], yy[r,c], z_bottom)
                v10 = (xx[r,c+1], yy[r,c+1], z_bottom)
                v01 = (xx[r+1,c], yy[r+1,c], z_bottom)
                v11 = (xx[r+1,c+1], yy[r+1,c+1], z_bottom)
                write_tri((0,0,-1), v00, v01, v10)
                write_tri((0,0,-1), v10, v01, v11)

        # 側面壁
        print("  側面壁...")
        wall_count = 0
        for r in range(nrows - 1):
            for c in range(ncols - 1):
                if not cell_mask[r, c]:
                    continue
                if r == 0 or not cell_mask[r-1, c]:
                    tl = (xx[r,c], yy[r,c], zz[r,c])
                    tr = (xx[r,c+1], yy[r,c+1], zz[r,c+1])
                    bl = (xx[r,c], yy[r,c], z_bottom)
                    br = (xx[r,c+1], yy[r,c+1], z_bottom)
                    n = calc_normal(bl, tl, br)
                    write_tri(n, bl, tl, br)
                    write_tri(n, br, tl, tr)
                    wall_count += 2
                if r == nrows - 2 or not cell_mask[r+1, c]:
                    tl = (xx[r+1,c], yy[r+1,c], zz[r+1,c])
                    tr = (xx[r+1,c+1], yy[r+1,c+1], zz[r+1,c+1])
                    bl = (xx[r+1,c], yy[r+1,c], z_bottom)
                    br = (xx[r+1,c+1], yy[r+1,c+1], z_bottom)
                    n = calc_normal(bl, br, tl)
                    write_tri(n, bl, br, tl)
                    write_tri(n, br, tr, tl)
                    wall_count += 2
                if c == 0 or not cell_mask[r, c-1]:
                    tu = (xx[r,c], yy[r,c], zz[r,c])
                    td = (xx[r+1,c], yy[r+1,c], zz[r+1,c])
                    bu = (xx[r,c], yy[r,c], z_bottom)
                    bd = (xx[r+1,c], yy[r+1,c], z_bottom)
                    n = calc_normal(bu, bd, tu)
                    write_tri(n, bu, bd, tu)
                    write_tri(n, bd, td, tu)
                    wall_count += 2
                if c == ncols - 2 or not cell_mask[r, c+1]:
                    tu = (xx[r,c+1], yy[r,c+1], zz[r,c+1])
                    td = (xx[r+1,c+1], yy[r+1,c+1], zz[r+1,c+1])
                    bu = (xx[r,c+1], yy[r,c+1], z_bottom)
                    bd = (xx[r+1,c+1], yy[r+1,c+1], z_bottom)
                    n = calc_normal(bu, tu, bd)
                    write_tri(n, bu, tu, bd)
                    write_tri(n, bd, tu, td)
                    wall_count += 2
        print(f"  側面壁三角形: {wall_count:,}")

        f.seek(count_pos)
        f.write(struct.pack('<I', counter[0]))

    _print_result(output_path, counter[0], model_width_mm, model_height_mm,
                  elev_range_mm + base_height_mm)


def grid_to_stl_rect(grid, lat_min, lat_max, lon_min, lon_max,
                     output_path, boundary_mask=None,
                     model_width_mm=160.0, z_exaggeration=2.0,
                     base_height_mm=3.0, smooth_boundary=True):
    """四角い台座付きSTLを生成（境界スムージング付き）"""
    nrows, ncols = grid.shape
    print(f"\nSTL生成中（四角い台座）... ({nrows} x {ncols})")

    mid_lat = (lat_min + lat_max) / 2.0
    lat_span_m = (lat_max - lat_min) * 111320.0
    lon_span_m = (lon_max - lon_min) * 111320.0 * np.cos(np.radians(mid_lat))

    scale = model_width_mm / lon_span_m
    model_height_mm = lat_span_m * scale

    if boundary_mask is not None and smooth_boundary:
        min_elev = max(0.0, float(np.min(grid[boundary_mask]))) if boundary_mask.any() else 0.0
        grid[~boundary_mask] = min_elev
        # Smoothstep境界ブレンド
        smooth_width = 5
        dist = np.zeros_like(grid, dtype=np.float32)
        dist[boundary_mask] = float(smooth_width)
        eroded = boundary_mask.copy()
        for d in range(1, smooth_width + 1):
            new_eroded = (np.roll(eroded, 1, axis=0) & np.roll(eroded, -1, axis=0) &
                          np.roll(eroded, 1, axis=1) & np.roll(eroded, -1, axis=1))
            edge = eroded & ~new_eroded
            dist[edge & (dist >= d)] = d
            eroded = new_eroded
        blend = np.clip(dist / smooth_width, 0.0, 1.0)
        blend = blend * blend * (3 - 2 * blend)
        grid = grid * blend + min_elev * (1 - blend)
        grid[~boundary_mask] = min_elev

    min_elev = float(np.min(grid))
    max_elev = float(np.max(grid))
    elev_range_mm = (max_elev - min_elev) * scale * z_exaggeration

    print(f"実寸: {lon_span_m/1000:.1f}km x {lat_span_m/1000:.1f}km")
    print(f"モデル: {model_width_mm:.1f} x {model_height_mm:.1f} x "
          f"{elev_range_mm + base_height_mm:.1f} mm")
    print(f"標高: {min_elev:.0f}m ～ {max_elev:.0f}m")

    x = np.linspace(0, model_width_mm, ncols, dtype=np.float32)
    y = np.linspace(model_height_mm, 0, nrows, dtype=np.float32)
    xx, yy = np.meshgrid(x, y)
    zz = ((grid - min_elev) * scale * z_exaggeration + base_height_mm).astype(np.float32)
    z_bottom = np.float32(0.0)

    total_tris = (nrows-1)*(ncols-1)*2*2 + 2*(2*(ncols-1) + 2*(nrows-1))
    print(f"三角形数: ~{total_tris:,}")

    with open(output_path, 'wb') as f:
        header = b'FGD DEM to STL - Rect Mode' + b'\0' * (80 - 26)
        f.write(header)
        count_pos = f.tell()
        f.write(struct.pack('<I', 0))

        write_tri, calc_normal, counter = _stl_helpers(f)

        # 上面
        for r in range(nrows - 1):
            if r % 300 == 0:
                print(f"  上面: {r}/{nrows-1}")
            for c in range(ncols - 1):
                v00 = (xx[r,c], yy[r,c], zz[r,c])
                v10 = (xx[r,c+1], yy[r,c+1], zz[r,c+1])
                v01 = (xx[r+1,c], yy[r+1,c], zz[r+1,c])
                v11 = (xx[r+1,c+1], yy[r+1,c+1], zz[r+1,c+1])
                n1 = calc_normal(v00, v10, v01)
                write_tri(n1, v00, v10, v01)
                n2 = calc_normal(v10, v11, v01)
                write_tri(n2, v10, v11, v01)

        # 底面
        print("  底面...")
        for r in range(nrows - 1):
            for c in range(ncols - 1):
                v00 = (xx[r,c], yy[r,c], z_bottom)
                v10 = (xx[r,c+1], yy[r,c+1], z_bottom)
                v01 = (xx[r+1,c], yy[r+1,c], z_bottom)
                v11 = (xx[r+1,c+1], yy[r+1,c+1], z_bottom)
                write_tri((0,0,-1), v00, v01, v10)
                write_tri((0,0,-1), v10, v01, v11)

        # 側面（4辺）
        print("  側面...")
        for c in range(ncols-1):
            r = nrows-1
            write_tri((0,-1,0), (xx[r,c],yy[r,c],z_bottom), (xx[r,c+1],yy[r,c+1],z_bottom),
                      (xx[r,c],yy[r,c],zz[r,c]))
            write_tri((0,-1,0), (xx[r,c+1],yy[r,c+1],z_bottom), (xx[r,c+1],yy[r,c+1],zz[r,c+1]),
                      (xx[r,c],yy[r,c],zz[r,c]))
        for c in range(ncols-1):
            r = 0
            write_tri((0,1,0), (xx[r,c],yy[r,c],z_bottom), (xx[r,c],yy[r,c],zz[r,c]),
                      (xx[r,c+1],yy[r,c+1],z_bottom))
            write_tri((0,1,0), (xx[r,c+1],yy[r,c+1],z_bottom), (xx[r,c],yy[r,c],zz[r,c]),
                      (xx[r,c+1],yy[r,c+1],zz[r,c+1]))
        for r in range(nrows-1):
            c = 0
            write_tri((-1,0,0), (xx[r,c],yy[r,c],z_bottom), (xx[r,c],yy[r,c],zz[r,c]),
                      (xx[r+1,c],yy[r+1,c],z_bottom))
            write_tri((-1,0,0), (xx[r+1,c],yy[r+1,c],z_bottom), (xx[r,c],yy[r,c],zz[r,c]),
                      (xx[r+1,c],yy[r+1,c],zz[r+1,c]))
        for r in range(nrows-1):
            c = ncols-1
            write_tri((1,0,0), (xx[r,c],yy[r,c],z_bottom), (xx[r+1,c],yy[r+1,c],z_bottom),
                      (xx[r,c],yy[r,c],zz[r,c]))
            write_tri((1,0,0), (xx[r+1,c],yy[r+1,c],z_bottom), (xx[r+1,c],yy[r+1,c],zz[r+1,c]),
                      (xx[r,c],yy[r,c],zz[r,c]))

        f.seek(count_pos)
        f.write(struct.pack('<I', counter[0]))

    _print_result(output_path, counter[0], model_width_mm, model_height_mm,
                  elev_range_mm + base_height_mm)


def _print_result(output_path, tri_count, w, h, z):
    file_size = os.path.getsize(output_path)
    print(f"\n{'='*50}")
    print(f"  完了！")
    print(f"  出力: {output_path}")
    print(f"  サイズ: {file_size / 1024 / 1024:.1f} MB")
    print(f"  三角形数: {tri_count:,}")
    print(f"  モデル: {w:.1f} x {h:.1f} x {z:.1f} mm")
    print(f"{'='*50}")


# =====================================================================
#  メイン
# =====================================================================

def main():
    parser = argparse.ArgumentParser(
        description='基盤地図情報 DEM → 3Dプリント用STL変換',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  # 糸魚川市の形だけのSTL（余白なし）
  python fgd_dem_to_stl.py ./dem_data -g N03.geojson -n 糸魚川 -o itoigawa.stl

  # 富士山周辺を四角い台座付きで
  python fgd_dem_to_stl.py ./dem_data -g N03.geojson -n 富士吉田 --rect -o fuji.stl

  # 境界なしで全データを変換
  python fgd_dem_to_stl.py ./dem_data -o all_terrain.stl

  # 高解像度版（解像度↑、ファイルサイズ↑）
  python fgd_dem_to_stl.py ./dem_data -g N03.geojson -n 糸魚川 --max-points 10000000
        """)
    parser.add_argument('input', help='DEMデータのフォルダパスまたはZIPファイル')
    parser.add_argument('-g', '--geojson', help='行政境界GeoJSONファイル')
    parser.add_argument('-n', '--name', help='切り抜く地域名（GeoJSON内で部分一致検索）')
    parser.add_argument('-o', '--output', default='terrain.stl', help='出力STLファイル名')
    parser.add_argument('-w', '--width', type=float, default=160.0,
                        help='モデルの幅(mm) [default: 160]')
    parser.add_argument('-z', '--z-exaggeration', type=float, default=2.0,
                        help='標高の誇張倍率 [default: 2.0]')
    parser.add_argument('-b', '--base-height', type=float, default=3.0,
                        help='底面の厚さ(mm) [default: 3.0]')
    parser.add_argument('--max-points', type=int, default=2_500_000,
                        help='グリッド最大ポイント数 [default: 2500000]')
    parser.add_argument('--rect', action='store_true',
                        help='四角い台座付きモードで出力')
    args = parser.parse_args()

    print("=" * 60)
    print("  基盤地図情報 DEM → 3Dプリント用STL")
    print("=" * 60)

    # 境界データ読み込み
    polygons = None
    bbox = None
    if args.geojson and args.name:
        print("\n─── 境界データ ───")
        polygons, bbox = load_boundary(args.geojson, args.name)
    elif args.geojson and not args.name:
        print("警告: --geojson を指定する場合は --name も必要です。境界なしで処理します。")
    elif args.name and not args.geojson:
        print("警告: --name を指定する場合は --geojson も必要です。境界なしで処理します。")

    # DEMファイル収集
    print("\n─── DEMデータ ───")
    xml_files = extract_all_xmls(args.input)
    if not xml_files:
        if os.path.isdir(args.input):
            xml_files = glob.glob(os.path.join(args.input, '**', '*.xml'), recursive=True)
    print(f"XMLファイル数: {len(xml_files)}")
    if not xml_files:
        print("エラー: XMLファイルが見つかりません")
        sys.exit(1)

    # グリッド構築
    print("\n─── グリッド構築 ───")
    grid, lat_min, lat_max, lon_min, lon_max, grid_lats, grid_lons = build_grid(
        xml_files, bbox, max_points=args.max_points
    )

    # STL生成
    print("\n─── STL生成 ───")
    if polygons and not args.rect:
        # 形状のみ
        boundary_mask = create_boundary_mask(grid_lats, grid_lons, polygons)
        grid_to_stl_shape(grid, lat_min, lat_max, lon_min, lon_max,
                          boundary_mask, args.output,
                          model_width_mm=args.width,
                          z_exaggeration=args.z_exaggeration,
                          base_height_mm=args.base_height)
    elif polygons and args.rect:
        # 四角い台座＋境界スムージング
        boundary_mask = create_boundary_mask(grid_lats, grid_lons, polygons)
        grid_to_stl_rect(grid, lat_min, lat_max, lon_min, lon_max,
                         args.output, boundary_mask=boundary_mask,
                         model_width_mm=args.width,
                         z_exaggeration=args.z_exaggeration,
                         base_height_mm=args.base_height)
    else:
        # 境界なし四角
        grid_to_stl_rect(grid, lat_min, lat_max, lon_min, lon_max,
                         args.output,
                         model_width_mm=args.width,
                         z_exaggeration=args.z_exaggeration,
                         base_height_mm=args.base_height)

    print(f"\n完成！")


if __name__ == '__main__':
    main()
