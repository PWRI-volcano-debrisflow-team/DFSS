---
affiliation: 土木研究所　火山・土石流チーム
author: 清水武志
date: 2023/1/19
---


# はじめに

古恵川のデータ（`../03_furuekawa` のこと。以下、オリジナルコードという。）を参考に、計算に必要な入出力ファイルをわかりやすくするためにフォルダ構成などを再整理したデータです。

以下にオリジナルコードとの相違点を示します。
ただし、前処理である地形データのGIS処理および後処理である可視化に関するプログラムは対象外としました。


## オリジナルコードの変更点について

`./bin` に fortran のソースコードを収める。オリジナルコードとそのパッチファイルを収めました。

下流域地形の前処理プログラムである `01_mk_watershedConfigurationFiles.f90` および `02_mk_floodplainConfigurationFiles.f90` は、
上流域地形の前処理プログラムと区別しやすくするために、ファイル名冒頭の数値を10番台に変更しました。

RR,DR,DFの変更点は、次の通りです。変更点の詳細は 各 `.patch` ファイルを参照してほしい。

+ 起動直後に `openMP` で使用コア数の標準出力を追加
+ このフォルダ構成で実行できるように入出力ファイルのパス名の変更


オリジナルコードから、以下のようにすれば、変更後のファイルが得られます。
```
patch -u < RR.path
```
これにより、`RR_ver1.0.f90` の内容が `RR_ver1.01.f90` と同様に変更されます。
元に戻すときは、
```
patch -R < RR.path
```
とします。

ソースコードのコンパイルを、例えば、`gfortran` で行う場合、以下のコマンドで実行できます。
```
gfortran 01_mk_topographyFiles.f90 -o 01_mk_topographyFiles.exe
gfortran 02_mk_watershedConfigurationFiles.f90 -o 02_mk_watershedConfigurationFiles.exe
gfortran 03_mk_streamConfigurationFiles.f90 -o 03_mk_streamConfigurationFiles.exe
gfortran 11_mk_floodplainConfigurationFiles.f90 -o 11_mk_floodplainConfigurationFiles.exe
gfortran 12_mk_streamFloodplainConnectionFiles.f90 -o 12_mk_streamFloodplainConnectionFiles.exe
gfortran -fopenmp RR_ver_1.01.f90 -o RR_ver_1.01_omp.exe
gfortran -fopenmp DR_ver_1.01.f90 -o DR_ver_1.01_omp.exe
gfortran -fopenmp DF_ver_1.01.f90 -o DF_ver_1.01_omp.exe
```


## サンプルデータについて

`./*.asc, ./*.txt` が DFSS の入力ファイルである。

`./input` フォルダはこのサンプルを動かすためには使用しないが、上の `./*.asc` 等を入力ファイルを作成するために、
ユーザが外部から収集・準備する以下のデータの保存を想定したフォルダである。

+ 流出解析を行う流域の下流末端の座標 ( `./input/coordinate.txt or ./input/outletPoint_furue.kml` ）
+ 国土地理院数値標高データ（10m B） `FG-GML-XXXX-XX-DEM10B.zip or .xml`
+ XRAINによる雨量データ
+ 下流域の計算範囲指定ファイル ( `./input/area_ll.kml` ）

`./input` フォルダにおける地形データと雨量データは利用規約に従って保存していません。

以下の処理を行なうことで、DFSS の計算に必要な入力ファイルが得られます。

+ 地形データは、[Open Source GRASS GIS](https://grass.osgeo.org/) を使ったGIS処理による測地系変換やリサンプリング、傾斜処理を実施して、`./dem_10m.asc, ./dir_10m.asc, ./dem_30m.asc, ./dir_30m.asc` に変換
+ 雨量データは、DFSS形式の雨量データに `./rain_Aso.txt` に変換

詳細は `../03_furuegawa/01_SET_BASIN` および `../03_furuegawa/02_MK_TOPO` を参照してください。なお、傾斜量のラスターデータ (`./dir_*m.asc`) は `GRASS GIS` による `r.watershed, r.fill.dir` によるアルゴリズムで得ることを想定しています。他の GIS ソフトウェアの処理では上手く動作しない例が確認されています。

さらに、計算にはパラメータファイル（`./topographyConfiguration.txt, ./RR_input.txt` 等）が必要です。



## 計算結果保存フォルダについて

`./bin` のコードは実行時に必要なフォルダの有無を判定して自動生成を行ないません。

そのため、計算前に `./output_RR, ./output_RR_misc, ./output_DR, ./output_DF` をあらかじめ作成しなければなりません。



## 並列計算で使用するコア数の設定と実行方法について

`openMP` による並列計算で使用する`CPU`のコア数は、windows のコマンドプロンプト（`cmd.exe`）や macOS や Linux のターミナル (`zsh`） によって、計算を実行する前に環境変数として設定する必要があります。

`./` がカレントディレクトリの場合に、コアを8つ使用するには以下のコマンドで実行できます。
```
set OMP_NUM_THREADS=8
bin\01_mk_topographyFiles.exe
bin\02_mk_watershedConfigurationFiles.exe
bin\03_mk_streamConfigurationFiles.exe
bin\11_mk_floodplainConfigurationFiles.exe
bin\12_mk_streamFloodplainConnectionFiles.exe
bin\RR_ver_1.01_omp.exe
bin\DR_ver_1.01_omp.exe
bin\DF_ver_1.01_omp.exe
```

なお、上に示した `cmd.exe` と `zsh` ではパス区切り文字が異なるため、使用環境にあわせて適宜変更してください。
