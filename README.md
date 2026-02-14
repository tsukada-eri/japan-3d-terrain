# japan-3d-terrain / 日本3D地形図メーカー

[日本語](#日本語) | [English](#english)

---

## 日本語

国土地理院の基盤地図情報DEMデータから、3Dプリント用STLファイルをワンコマンドで生成するツールです。

QGISなどのGISソフトウェアは不要。Python + NumPy だけで動きます。

### 特徴

- DEM5A/5B/5C（5mメッシュ）とDEM10B（10mメッシュ）の混在データを自動統合
- GeoJSON行政境界で任意の市区町村の形に切り抜き
- 市域形状のみのSTL（余白なし）or 四角い台座付きを選択可能
- ZIPファイル（入れ子含む）の自動展開
- メモリ制限付き自動ダウンサンプル

### クイックスタート

```bash
pip install numpy

# 糸魚川市の形で切り抜き
python fgd_dem_to_stl.py ./dem_data \
  -g N03-23_15_230101.geojson \
  -n 糸魚川 \
  -o itoigawa.stl \
  -w 160 -z 2.0

# 四角い台座付き
python fgd_dem_to_stl.py ./dem_data \
  -g boundary.geojson -n 糸魚川 --rect -o itoigawa_rect.stl

# 境界なし（全データ）
python fgd_dem_to_stl.py ./dem_data -o terrain.stl
```

### オプション

| オプション | 説明 | デフォルト |
|-----------|------|----------|
| `input` | DEMデータのフォルダまたはZIP | （必須） |
| `-g` | 行政境界GeoJSON | なし |
| `-n` | 切り抜く地域名（部分一致） | なし |
| `-o` | 出力STLファイル名 | `terrain.stl` |
| `-w` | モデル幅 (mm) | `160` |
| `-z` | 標高の誇張倍率 | `2.0` |
| `-b` | 底面の厚さ (mm) | `3.0` |
| `--max-points` | グリッド最大ポイント数 | `2500000` |
| `--rect` | 四角い台座モード | off |

### 必要なデータの入手先

| データ | 入手先 | 用途 |
|--------|--------|------|
| 基盤地図情報 DEM | https://fgd.gsi.go.jp/download/menu.php | 標高データ（無料・要登録） |
| 行政区域データ (N03) | https://nlftp.mlit.go.jp/ksj/ | 市区町村の境界線（無料） |

詳しい手順は [docs/manual_ja.md](docs/manual_ja.md) を参照してください。

### AIアシスタントとの連携

LLM（大規模言語モデル）にこのスクリプトを渡して作業を依頼できます。
プロンプトテンプレートは [docs/ai_prompt.md](docs/ai_prompt.md) を参照してください。

### ライセンス

MIT License - 自由に利用・改変・再配布できます。

---

## English

Generate 3D-printable STL terrain models from Japan's Geospatial Information Authority (GSI) DEM data with a single command.

No GIS software required. Just Python + NumPy.

### Features

- Automatic merging of DEM5A/5B/5C (5m mesh) and DEM10B (10m mesh) with bilinear resampling
- Clip to any municipality shape using GeoJSON administrative boundaries
- Shape-only STL (no bounding box) or rectangular base plate mode
- Auto-extraction of nested ZIP files
- Memory-safe auto-downsampling

### Quick Start

```bash
pip install numpy

# Clip to Itoigawa city shape
python fgd_dem_to_stl.py ./dem_data \
  -g N03-23_15_230101.geojson \
  -n 糸魚川 \
  -o itoigawa.stl \
  -w 160 -z 2.0

# Rectangular base plate mode
python fgd_dem_to_stl.py ./dem_data \
  -g boundary.geojson -n 糸魚川 --rect -o itoigawa_rect.stl

# No boundary (full extent)
python fgd_dem_to_stl.py ./dem_data -o terrain.stl
```

### Options

| Option | Description | Default |
|--------|-------------|---------|
| `input` | DEM data folder or ZIP file | (required) |
| `-g` | Administrative boundary GeoJSON | none |
| `-n` | Area name to clip (partial match) | none |
| `-o` | Output STL filename | `terrain.stl` |
| `-w` | Model width (mm) | `160` |
| `-z` | Elevation exaggeration | `2.0` |
| `-b` | Base height (mm) | `3.0` |
| `--max-points` | Max grid points | `2500000` |
| `--rect` | Rectangular base mode | off |

### Required Data Sources

| Data | Source | Purpose |
|------|--------|---------|
| FGD DEM | https://fgd.gsi.go.jp/download/menu.php | Elevation data (free, registration required) |
| Admin boundaries (N03) | https://nlftp.mlit.go.jp/ksj/ | Municipality borders (free) |

See [docs/manual_ja.md](docs/manual_ja.md) for detailed instructions (Japanese).

### Using with AI Assistants

You can hand this script to an LLM and ask it to run the conversion for you.
See [docs/ai_prompt.md](docs/ai_prompt.md) for a ready-made prompt template.

### License

MIT License
