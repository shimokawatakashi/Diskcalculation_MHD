# diskcalculation

原始惑星系円盤の長期進化を、ガス・ダスト・磁場・温度を同時に解いて計算する 1 次元 (動径方向) コード。
鉛直方向の構造 (電離度・磁気拡散係数) は各半径で鉛直格子を用いて求める。

## 研究の目的

原始惑星系円盤で観測されるリング・ギャップ構造が、磁場の輸送と非理想 MHD 効果 (オーム散逸・両極性拡散) によってどのように形成されるかを明らかにする。
そのために、次の過程を同時に解く円盤の大局的な進化計算を行う。

- ダストを含む電離平衡から求めた磁気拡散係数による、鉛直磁場 (磁束) の動径方向の輸送
- 磁場と磁気拡散で決まる MRI 乱流の強さ (α_MRI) と、円盤風による角運動量・質量の損失 (α_DW, Σ̇_w)
- ダストの合体成長・破砕・動径方向の移動と、それが電離度 (磁気拡散) に与える影響
- 粘性加熱・ジュール加熱・円盤風加熱と輻射冷却で決まる温度

宇宙線電離率や初期プラズマ β などの条件を変えたモデルを比べ、擬似観測 (ダスト連続波の強度分布) を作って観測と比較する。

## コンパイル

`g++` (C++17) と OpenMP を使う。リポジトリ直下で

```sh
make            # bin/run を作成 (src/*.cpp と Marchand2022/src の 2 ファイルをコンパイル)
make clean      # obj/*.o と bin/run を削除
```

コンパイラやフラグは `Makefile` の先頭 (`CXX`, `CXXFLAGS`, `LDFLAGS`) で変更する。MPI を使うときは `CXX = mpic++` と、コメントアウトしてあるリンク規則に切り替える。

## 実行

入力ファイルは `input_folder/` に置く (`input_folder/input` が全パラメータの説明付きの例)。

```sh
make run ARGS=input_folder/input                  # 新規計算 (OMP_NUM_THREADS=4)
./bin/run input_folder/input                      # 直接実行する場合
```

途中の出力から再開する場合:

```sh
make restart ARGS=input_folder/input BASE=/mnt/hdd/shimokawa RUN=test YEAR=320000
# = ./bin/run -r <入力> <BASE/RUN/outYEAR> <BASE/RUN/outYEAR_gas> <BASE/RUN_B/outYEAR_B> <YEAR>
```

長時間の計算はスクリプトで自動的に区切って再開できる。

```sh
make supervised ARGS=input_folder/input   # scripts/run_supervised.sh : SEGMENT_HOURS ごとに区切り, 最新の出力から自動で再開
make daytime    ARGS=input_folder/input   # scripts/run_daytime.sh    : 毎日 WINDOW_START–WINDOW_END の間だけ実行
```

主な入力パラメータ:

| キー | 意味 |
|---|---|
| `outdir` | 出力先。磁場は `outdir` に `_B` を付けたディレクトリに出力 |
| `diskin`, `diskout`, `nrdisk`, `radial_grid_type` | 動径格子 (内側・外側半径 [au], 格子数, 格子の種類) |
| `z_max`, `nzdisk` | 鉛直格子 (上端 [H_g], 格子数) |
| `sigmag0_`, `cut_off`, `star_mass` | 初期ガス面密度の規格化, カットオフ半径, 中心星質量 |
| `amin`, `amax`, `ndust`, `nd_lin`, `rhodi`, `fdg` | ダストのサイズビンと内部密度, ダスト・ガス比 (0 でダストなし) |
| `beta_ini`, `B_out`, `eta_max`, `zeta_CR0` | 初期プラズマ β, 磁場の下限, η の上限, 宇宙線電離率 |
| `alpha_t`, `alpha_DW`, `wind_alpha_model`, `lambda_w`, `Cw` | 乱流・円盤風のパラメータ |
| `growth`, `fragmentation`, `thermal` | 合体成長・破砕・温度計算の on/off |
| `use_implicit_mag`, `mag_theta`, `use_implicit_coag`, `use_full_implicit_coag`, `face_flux_scheme` | 数値解法の選択 |
| `cfl`, `tend` | CFL 数, 終了時刻 [yr] |

計算の種類は `growth` と `fragmentation` で決まる (`DustDisk::CalculationRun`)。

| `growth` | `fragmentation` | 呼ばれる関数 |
|---|---|---|
| true | true | `UpdateDustDiskFragOn` (合体成長 + 破砕) |
| true | false | `UpdateDustDisk` (合体成長のみ。**工事中: エラーで停止する**) |
| false | – | `UpdateNoGrowthDustDisk` (成長なし。ダストは移流・拡散・円盤風による損失のみ) |

## 出力ファイル

`outdir/out<年>` (ダスト):
1–4 行目が `m_` (質量ビン), `mb_` (ビン境界), `a_` (半径ビン), `ab_` (ビン境界)、空行のあと各行が `r [au]` と各ビンの面密度 `sigmad_active [g/cm²]`。

`outdir/out<年>_gas` (列番号は 0 から):

| 列 | 量 | 列 | 量 |
|---|---|---|---|
| 0 | r [au] | 13 | 中心面のガス圧 P |
| 1 | Σ_g | 14 | 中心面の磁気圧 P/β |
| 2 | Σ_d (全ダスト) | 15 | dN/dr (粘性トルク) |
| 3 | Σ_vapor (昇華したダスト) | 16 | α_GI |
| 4 | v_g,r | 17 | α_φz (円盤風) |
| 5 | v_g,φ | 18 | α_tot |
| 6 | v_wind | 19 | η_P |
| 7 | v_d,r (最小ビン, 合計) | 20–24 | 円盤風モデルの z_dead, z_surf, Δ, Am_eff, α_DW |
| 8 | v_d,r (最小ビン, ドリフト) | 25–28 | (`face_flux_scheme = true` のみ) r_face, v_g,r,face, Δv_φ,face, v_Bdiff,face |
| 9 | v_d,r (最小ビン, 拡散) | | |
| 10 | ダスト・ガス比 (最小ビン) | | |
| 11 | dΣ_g/dt (円盤風) | | |
| 12 | 中心面温度 T | | |

`outdir_B/out<年>_B`: 0: r, 1: α_MRI, 2: B_z, 3: η_O, 4: η_A, 5: 中心面のプラズマ β。

## プロット

`plot/` の Jupyter ノートブックで解析する。

1. `plot/thesis_common.ipynb` の「解析設定」で、データの場所 (`BASE_DIR`)、モデル (`MODEL`)、表示範囲、図のラベル (`model_names`, `FIG_TAGS`) を指定する。
2. 各ノートブックを先頭のセルから実行する。先頭のセルで `thesis_common.ipynb` が読み込まれる。

| ノートブック | 内容 |
|---|---|
| `analysis_snapshot.ipynb` | 時刻を固定した動径分布 (面密度, 速度, ストークス数, タイムスケール, 磁束, β など) |
| `analysis_spacetime.ipynb` | (r, t) 図と時間進化 (磁場の輸送, 質量の時間進化, リング・ギャップ) |
| `analysis_obs_theory.ipynb` | 擬似観測によるリング・ギャップと, 解析解との比較 |
| `compare_beta_ini.ipynb` | 複数モデル (宇宙線電離率など) の比較。モデルと色は `RUN_COLORS` で指定 |
| `anim_profiles.ipynb` | Σ, T, B の動径分布のアニメーション (GIF) |
| `gap_detection.ipynb` | 強度分布からリング・ギャップの位置を求める (他のノートブックから `run_nb` で呼ぶ) |

擬似観測のダスト吸収係数は `optool/` の `optool -dsharp` で計算するので、事前に `optool` をビルドしておく。図は `./fig` (解析用) と `./phd_fig` (論文用) に保存される。

## ソースコードと主な関数

### `src/main.cpp`
入力ファイルを読み、`DustDisk` を初期化 (`Initialization` または再開時は `Initialization_Restart`) して `CalculationRun` を呼ぶ。

### `src/dust_disk.cpp` (`DustDisk` クラス)

**入力と初期化**
- `InputParameter` : 入力ファイルからパラメータを読む。出力ディレクトリを作り入力ファイルをコピーする
- `Initialization` / `Initialization_Restart` : 格子・初期分布の設定 / 出力ファイルからの再開
- `SetRadiusGrid`, `SetVerticalGrid` : 動径・鉛直格子
- `SetDiffrentialCoefficient_R`, `_R_Half`, `For2order_R`, `_Z`, `SetZeroDifferentialCoeffient_R` : 不等間隔格子の差分係数
- `SetBoundaryConditionForDisk`, `SetBoundaryConditionForMagneticfield` : 境界条件
- `SetDustGrid` : ダストの質量・半径ビン
- `SetDiskGasModel` : 初期のガス円盤 (面密度, 温度, 磁場)
- `SetInitialDustDistribution`, `CalcSizeForTargetStokes`, `CalcInitialStokesLimit` : 初期ダスト分布 (monomer / 冪分布)
- `CalcXrayIonizationrate`, `CalcCosmicRayIonizationrate` : X 線・宇宙線による電離率
- `Iteration` : 温度・磁気拡散・乱流を反復して初期分布を計算

**物理量の更新**
- `UpdatePhysicalParameter` : 音速, スケールハイト, 密度, 圧力など
- `DetermineVaporDust` : スノーラインの内側でのダストの昇華・凝結
- `CalcStokesNumberAndDustScaleHeight` : ストークス数とダストのスケールハイト
- `UpdateResistivity` : 各 (r, z) で電離平衡を解き η_O, η_A を求める
- `UpdateTurbulent` : α_tot, 粘性, ダストの拡散係数
- `CalcWindAlpha`, `CalcWindTorque_Suzuki2016`, `CalcWindDrivenAlpha_Enomoto` : 円盤風による α_DW (`wind_alpha_model` で切り替え)
- `CalcWindMassLoss_Tabone2022`, `CalcWindMassLoss_Riols2019` : 円盤風による質量損失率
- `CalcMagneticFieldPressure`, `CalculateViscousTorque`, `CalculateDust_to_Gas_massratio`, `CalculateDiskCoeffcient` : 磁気圧, 粘性トルク, ダスト・ガス比, 速度の式の係数
- `CalcGasVelocity_TM2018`, `CalcDustVelocity_TM2018`, `CalcDustVelocity` : ガス・ダストの動径・方位角方向の速度
- `CalcGasFlux`, `CalcDustFlux` : 面密度のフラックス
- `ApplyDustDiffusionImplicit` : ダストの乱流拡散 (陰解法)

**時間発展と出力**
- `CalculationRun` : `growth` / `fragmentation` に応じて時間発展の関数を選ぶ
- `UpdateDustDiskFragOn` : 合体成長 + 破砕ありの時間発展 (メインの計算)
- `UpdateDustDisk` : 合体成長のみの時間発展 (工事中)
- `UpdateNoGrowthDustDisk` : ダスト成長なしの時間発展
- `OutPut` : ダスト・ガス・磁場のファイルを書き出す

### `src/dust_growth.cpp` (`DustDisk` クラス, ダストの合体成長・破砕)
- `CalcCoagulationAndFragmentationArrays` : 衝突後の質量の配分 (ビン番号と割合) の表
- `CalcDustRelativeVelocity` : 相対速度 (ブラウン運動, 乱流, 動径・方位角・鉛直ドリフト)
- `CoagulateChangeRate` : 合体のみの変化率 (旧版, 現在は未使用)
- `CoagulateChangeRateFragmentation` : 合体・破砕による変化率 (陽解法)
- `CoagulateChangeRateFragmentation_Implicit`, `_Implicit_Cell` : 損失項を陰的に扱う変化率 (半陰解法)
- `CoagulateChangeRateFragmentation_FullImplicit_Cell` : Jacobian を含む変化率 (完全陰解法)
- `ApplyCoagulationAndCorrectMass`, `_Implicit`, `_FullImplicit` : 各解法で面密度を更新し, 質量保存を補正する
- `EstimateDtCoagulation`, `EstimateDtGain` : 陽解法・半陰解法の時間刻みの制限

### `src/face_scheme.cpp` (`DustDisk` クラス, `face_flux_scheme = true` のとき)
セル境界の速度を使い, MUSCL で移流を解く。
`SetupFaceScheme`, `CalcGasFaceVelocity`, `AdvanceGasFaceScheme`, `CalcDustFaceVelocity`, `CalcDustFluxFace`, `CalcVaporFluxFace`, `AdvanceDustAdvectionFace`, `AdvanceVaporFace`, `GasAdvectionDtFace`, `DustAdvectionDtFace`。

### `src/alpha.cpp` (`CalcAlphaParameter` クラス, 電離平衡と磁気拡散係数)
- `SetParameter`, `SetParameterForAlpha` : 電離率, 温度, ガス・ダストの数密度, 磁場などを設定
- `GetChargestate` : ダスト表面電位 φ についての電荷中性条件を Newton 法で解く (`total_dustcharge`, `total_ion_numberdensity`, `total_electron_numberdensity` と各 `diff_*` が補助関数)
- `CalcThermalIonization` : カリウムの熱電離 (Saha の式, T > 1000 K)
- `CalcResistivity`, `CalcResistivity_ForThermal` : η_O, η_A (後者は熱電離を含む)
- `GetAlphaParameter` : 磁気拡散係数から α を求める

### `src/magnetic_field1D.cpp` (`Magnetic_Field2D` クラス, 鉛直磁場の輸送)
- `SetRadiusGridParameter`, `SetVerticalGridParameter`, `Set_Differential_Coefficient`, `SetDiskParameter`, `SetMagneticFieldParameter` : 格子と円盤・磁場の量を設定
- `CalculateDeltaTime` : 磁場の陽解法の時間刻み
- `Calc_Mag`, `Calc_Mag_Implicit` : 誘導方程式の時間発展 (移流 `Calculate_Advection_B` + 拡散 `Calculate_Diffuse_B`, 陰解法は θ 法の三重対角)
- `CalculationMRIturbulence : β, エルサッサー数から α_MRI を求める 

### `src/energy.cpp` (`Energy_Calc` クラス, 温度)
- `SetGridParameter`, `SetDiskParameter`, `SetMagneticField` : 格子と円盤・磁場の量を設定
- `RosselandOpacity` : Rosseland 平均吸収係数
- `CalcViscousHeating`, `CalcJouleheating`, `CalcWindHeating` : 粘性加熱, ジュール加熱, 円盤風加熱
- `CalcTemperature`, `GetTemperature`, `GetTemperature_Initialization`, `GetTemperature_TimeEvolution` : 加熱と冷却 (照射を含む) から温度を求める
- `Calculate_Dt` : 温度計算の時間刻み

### `src/utils.cpp`
`InputConfigure` (入力ファイルの読み込み: `ReadFile`, `Get`, `Contains`), `Split`, `FileCopy`, `ConvertStrToBool`。

### その他
- `src/physical_constants.hpp`, `src/typedef.hpp` : 物理定数と型の定義
- `Marchand2022/` : 磁気拡散係数の計算 (Marchand et al. 2022 の方法)。`calc_resistivity.cpp` と `calculation_tool.cpp` を本体と一緒にコンパイルする
- `optool/` : ダストの吸収係数の計算 (プロット用)
- `scripts/` : 長時間計算の自動再開スクリプト

## 参考文献 (コード中で引用しているもの)

**ダストの合体成長・破砕・運動**
- Brauer, F., Dullemond, C. P., & Henning, Th. 2008 — ダストのスケールハイト, 完全陰解法
- Birnstiel, T., et al. 2010 — ダストのスケールハイト, 鉛直ドリフト
- Birnstiel, T., Klahr, H., & Ercolano, B. 2012 — 破砕・ドリフトで決まる最大サイズ (初期分布)
- Garaud, P., et al. 2013 — 相対速度の分布を考えた合体・破砕の衝突率
- Kawasaki, Y., et al. 2022 — 破砕モデル, 質量の配分
- Okuzumi, S., et al. 2012 — 方位角方向のドリフトによる相対速度
- Ormel, C. W., & Cuzzi, J. N. 2007 — 乱流による相対速度
- Stammler, S. M., & Birnstiel, T. 2022 (DustPy) — 陰解法の時間刻みの採否

**電離・磁気拡散係数**
- Marchand, P., et al. 2016; Marchand, P. 2021; Marchand, P., et al. 2022 — 電離平衡と磁気拡散係数, 衝突率
- Okuzumi, S. 2009 — ダストを含む電離平衡
- Pinto, C., & Galli, D. 2008 — 荷電粒子と中性ガスの衝突率
- Sano, T., & Turner, N. J. 2008 — 電離率
- Tsukamoto, Y., & Okuzumi, S. 2022 — ダストを含む電離平衡 (電荷中性条件)
- Umebayashi, T., & Nakano, T. 1990 — 電子の付着確率
- Umebayashi, T., & Nakano, T. 2009 — 電離率

**磁場・乱流・円盤風**
- Delage, T. N., et al. 2022 — MRI 乱流の α とエルサッサー数による条件
- Enomoto, Mori, & Okuzumi 2026 (PASJ) — 円盤風による降着の自己無撞着な α_DW
- Riols, A., et al. 2019 — 円盤風による質量損失, 磁気粘性
- Sano, T., & Miyama, S. M. 1999 — MRI の最大成長波長
- Suzuki, T. K., et al. 2016 — 円盤風のトルクと質量損失
- Tabone, B., et al. 2022 — 円盤風による質量損失 (磁気レバーアーム λ)

**温度**
- Kunitomo, M., et al. 2020 — 円盤の温度の式と吸収係数

**プロット (擬似観測)**
- Aso, Y., & Machida, M. N. 2020; Yamaguchi, T., et al. 2024 — リング・ギャップの判定
- Birnstiel, T., et al. 2018 (DSHARP) — ダストの吸収係数 (`optool -dsharp`)

---

Created by Claude Code (2026-11-21)
