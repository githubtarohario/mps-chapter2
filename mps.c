/*=====================================================================
  mps.c
  (c) Kazuya SHIBATA, Kohei MUROTANI and Seiichi KOSHIZUKA (2014)

   Fluid Simulation Program Based on a Particle Method (the MPS method)
   Last update: May 21, 2014
=======================================================================*/
/*=====================================================================
 【このプログラムの概要】

   粒子法の一種である MPS 法 (Moving Particle Semi-implicit method /
   半陰的粒子法) によって非圧縮性流体の流れを解くプログラムである。

   計算対象は「ダムブレイク (ダム崩壊) 問題」。
   箱の左側に水柱 (幅 0.25m × 高さ 0.50m) を置き、重力で崩れて
   右方向へ流れ広がる様子を 2 秒間 (FINISH_TIME) 計算する。

 【MPS 法の 1 ステップの流れ (半陰的アルゴリズム)】

   (1) 陽的計算 (explicit part) … 圧力以外の外力を先に効かせる
        calGravity()      : 重力加速度を与える
        calViscosity()    : 粘性項 (ラプラシアンモデル) を加える
        moveParticle()    : 速度・位置を仮に進める (仮の速度 u*, 仮の位置 r*)
        collision()       : 粒子どうしがめり込んだ場合の反発処理

   (2) 陰的計算 (implicit part) … 非圧縮性を満たす圧力を求める
        calPressure()     : 圧力のポアソン方程式 [A]{P}={b} を解く
                            └ calNumberDensity()  : 粒子数密度 n を数える
                            └ setBoundaryCondition(): 自由表面/内部の判定
                            └ setSourceTerm()     : 右辺ベクトル {b}
                            └ setMatrix()         : 係数行列 [A]
                            └ solveSimultanious…  : ガウスの消去法で求解
                            └ removeNegativePressure(): 負圧をカット
                            └ setMinimumPressure(): 近傍最小圧力を記録

   (3) 圧力による修正 (correction part)
        calPressureGradient()          : 圧力勾配から加速度を計算
        moveParticleUsingPressureGradient() : 速度・位置を修正して確定

 【粒子の種類 (ParticleType)】
     FLUID      ( 0) : 流体粒子。運動方程式を解く対象。
     WALL       ( 2) : 壁粒子。動かないが圧力計算には参加する。
     DUMMY_WALL ( 3) : ダミー壁粒子。壁の外側に並べ、壁粒子の粒子数密度
                       が不足しないようにするためだけの粒子。
                       圧力は解かない (ノイマン境界の代用)。
     GHOST      (-1) : 計算対象外となった粒子 (本プログラムでは未使用)。

 【出力ファイル】
     output_%04d.prof   : テキスト形式。時刻・粒子数・各粒子の状態量。
     particle_%04d.vtu  : VTK 形式。ParaView 等で可視化できる。
     いずれも実行ファイルのカレントディレクトリに OUTPUT_INTERVAL
     ステップごとに書き出される。

 【Visual C++ (VC++) でのビルド方法】

   ● コマンドライン (開発者コマンドプロンプト) の場合
         cl /W3 /O2 mps.c /Fe:mps.exe
         mps.exe

   ● Visual Studio の IDE の場合
         1. [ファイル]-[新規作成]-[プロジェクト] で
            「空のプロジェクト」(C++) を作成する。
         2. 「ソースファイル」に本 mps.c を追加する。
            ※拡張子が .c なので C コンパイラとしてビルドされる。
         3. 構成を Release / x64 にして [ビルド]-[ソリューションのビルド]。
         4. 出力ファイルはカレントディレクトリ
            (IDE から実行した場合はプロジェクトフォルダ) に作られる。

   ※本ファイルは UTF-8 (BOM 付き) で保存してある。BOM が無いと
     VC++ が日本語コメントを CP932 と誤認して文字化けするため、
     編集して保存し直すときも UTF-8 (BOM 付き) を維持すること。

 【VC++ 向けに行った修正 (オリジナルからの変更点)】

   1) 先頭に _CRT_SECURE_NO_WARNINGS を定義した。
      VC++ では sprintf() / fopen() が「安全でない関数」とされ、
      警告 C4996 が出る (「/sdl」や「警告をエラーとして扱う」設定では
      ビルドが失敗する)。この定義でその警告を抑止している。

   2) fopen() の戻り値が NULL (ファイルを開けなかった) の場合の
      チェックを追加した。オリジナルはチェックが無いため、書き込み
      禁止フォルダで実行すると NULL ポインタ参照で異常終了する。

   3) main() の未使用引数 argc / argv を (void) キャストして、
      警告レベル /W4 での C4100 警告を抑止した。

   ※計算内容 (数式・アルゴリズム) は一切変更していない。

 【メモリについての注意】
   係数行列 CoefficientMatrix は ARRAY_SIZE×ARRAY_SIZE の静的配列で、
   5000×5000×8byte = 約 200MB を占有する。静的領域 (.bss) なので
   VC++ の x86/x64 いずれでもビルド・実行できるが、ARRAY_SIZE を
   大きくしすぎるとメモリ不足になるので注意すること。
=======================================================================*/

/* VC++ で sprintf / fopen を使うと警告 C4996 が出るため、
   それを抑止する。必ず #include より前に書くこと。            */
#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>      /* printf, fprintf, fopen, fclose, sprintf 用  */
#include <stdlib.h>     /* exit 用                                      */
#include <math.h>       /* sqrt 用                                      */
#include <string.h>     /* 文字列操作用 (オリジナルのまま残してある)    */

/*---------------------------------------------------------------------
  【計算条件の設定 (2 次元計算用)】
  ここを書き換えれば計算条件を変更できる。
---------------------------------------------------------------------*/
/* ★現在は「3 次元計算」に設定されている★
   2 次元に戻すには、下の 2 次元用 4 行の行頭の // を外し、
   さらにその下の 3 次元用 4 行の行頭に // を付ける。          */

/* ----- 2 次元計算用（現在は無効） ----- */
//#define DIM                  2       /* 空間次元数。2 なら 2 次元計算   */
//#define PARTICLE_DISTANCE    0.025   /* 初期粒子間距離 l0 [m]。
//                                        粒子を並べる間隔であり、
//                                        いわゆる「解像度」に相当する。
//                                        小さくすると精度は上がるが
//                                        粒子数が増えて計算時間が激増する */
//#define DT                   0.001   /* 時間刻み Δt [s]。
//                                        大きすぎると計算が発散する      */
//#define OUTPUT_INTERVAL      20      /* 何ステップごとにファイル出力するか
//                                        20 なら 20×DT = 0.02 秒ごと     */

/* for three-dimensional simulation */
/* ----- 3 次元計算用（現在有効） -----
   2 次元より粒子間距離を粗く (0.075) してある。
   3 次元では粒子数が一気に増えるため、2 次元と同じ 0.025 にすると
   粒子数が数十万個になり現実的な時間で計算できないからである。
   Δt を 0.003 と大きくしているのも同じ理由による。            */
#define DIM                  3       /* 空間次元数。3 なら 3 次元計算   */
#define PARTICLE_DISTANCE    0.075   /* 初期粒子間距離 l0 [m]           */
#define DT                   0.003   /* 時間刻み Δt [s]                 */
#define OUTPUT_INTERVAL      2       /* 2 なら 2×DT = 0.006 秒ごとに出力 */

/*---------------------------------------------------------------------
  【物理定数・計算パラメータ】
---------------------------------------------------------------------*/
#define ARRAY_SIZE           5000    /* 配列に確保する粒子数の上限。
                                        実際の粒子数がこれを超えると
                                        配列外アクセスとなるので注意   */
#define FINISH_TIME          2.0     /* 計算終了時刻 [s]                */
#define KINEMATIC_VISCOSITY  (1.0E-6)/* 動粘性係数 ν [m^2/s] (20℃の水) */
#define FLUID_DENSITY        1000.0  /* 流体密度 ρ [kg/m^3] (水)        */
#define G_X  0.0                     /* 重力加速度の x 成分 [m/s^2]     */
#define G_Y  -9.8                    /* 重力加速度の y 成分 [m/s^2]。
                                        下向き (−y 方向) が重力         */
#define G_Z  0.0                     /* 重力加速度の z 成分 [m/s^2]     */

/* 影響半径 re : 重み関数 w(r) が 0 でなくなる距離。
   この距離より近い粒子だけを「近傍粒子」として相互作用させる。
   物理量の種類ごとに最適な値が異なるため 3 種類定義している。        */
#define RADIUS_FOR_NUMBER_DENSITY  (2.1*PARTICLE_DISTANCE)
                                     /* 粒子数密度 n を計算するときの re */
#define RADIUS_FOR_GRADIENT        (2.1*PARTICLE_DISTANCE)
                                     /* 勾配モデル (圧力勾配) 用の re   */
#define RADIUS_FOR_LAPLACIAN       (3.1*PARTICLE_DISTANCE)
                                     /* ラプラシアンモデル (粘性・圧力
                                        ポアソン方程式) 用の re。
                                        勾配より広くとるのが定石       */

#define COLLISION_DISTANCE         (0.5*PARTICLE_DISTANCE)
                                     /* 衝突判定距離。粒子間距離がこれ
                                        より小さくなったら「めり込み」
                                        とみなし反発させる             */
#define THRESHOLD_RATIO_OF_NUMBER_DENSITY  0.97
                                     /* 自由表面判定のしきい値 β。
                                        n < β×n0 なら自由表面粒子。
                                        表面では近傍粒子が欠けるので
                                        粒子数密度が小さくなる         */
#define COEFFICIENT_OF_RESTITUTION 0.2
                                     /* 衝突の反発係数 e (0:完全非弾性,
                                        1:完全弾性)                    */
#define COMPRESSIBILITY (0.45E-9)    /* 水の圧縮率 [1/Pa]。
                                        圧力方程式の対角項に足して
                                        計算を安定化させる             */
#define EPS             (0.01 * PARTICLE_DISTANCE)
                                     /* 微小量。浮動小数点の丸め誤差で
                                        領域判定を誤らないようにする   */
#define ON              1            /* 論理値「真」                   */
#define OFF             0            /* 論理値「偽」                   */
#define RELAXATION_COEFFICIENT_FOR_PRESSURE 0.2
                                     /* 圧力ポアソン方程式の右辺に掛け
                                        る緩和係数 γ。1 未満にすると
                                        圧力振動が抑えられる           */

/*---------------------------------------------------------------------
  【粒子の種類 ParticleType[] に入る値】
---------------------------------------------------------------------*/
#define GHOST  -1                    /* 計算対象外の粒子               */
#define FLUID   0                    /* 流体粒子                       */
#define WALL    2                    /* 壁粒子 (圧力を解く)            */
#define DUMMY_WALL  3                /* ダミー壁粒子 (圧力を解かない)  */

/*---------------------------------------------------------------------
  【圧力計算での境界条件 BoundaryCondition[] に入る値】
---------------------------------------------------------------------*/
#define GHOST_OR_DUMMY  -1           /* 圧力を解かない粒子
                                        (ゴースト粒子・ダミー壁粒子)   */
#define SURFACE_PARTICLE 1           /* 自由表面粒子。
                                        圧力 P=0 のディリクレ境界      */
#define INNER_PARTICLE   0           /* 内部粒子。圧力を未知数として解く */

/*---------------------------------------------------------------------
  【ディリクレ境界の連結性チェック用のフラグ】
   圧力ポアソン方程式は、粒子のかたまりの中に P=0 を与える点
   (自由表面粒子) が 1 つも無いと解が一意に決まらない (特異行列)。
   そこで「自由表面粒子とつながっているか」を幅優先的に調べる。
---------------------------------------------------------------------*/
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0 /* まだつながっていない    */
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1 /* つながった (探索待ち)   */
#define DIRICHLET_BOUNDARY_IS_CHECKED       2 /* 探索済み                */

/*=====================================================================
  【関数プロトタイプ宣言】
  すべての関数は引数なし (void)、戻り値なし (void) で、
  グローバル変数を介してデータをやり取りする設計になっている。
  例外は weight() のみで、引数 2 つ・戻り値 double である。
=====================================================================*/
void initializeParticlePositionAndVelocity_for2dim( void ); /* 2次元の初期配置 */
void initializeParticlePositionAndVelocity_for3dim( void ); /* 3次元の初期配置 */
void calConstantParameter( void );        /* 定数パラメータの事前計算   */
void calNZeroAndLambda( void );           /* 基準粒子数密度 n0 と λ     */
double weight( double distance, double re ); /* 重み関数 w(r)           */
void mainLoopOfSimulation( void );        /* 時間発展のメインループ     */
void calGravity( void );                  /* 重力加速度の設定           */
void calViscosity( void );                /* 粘性項の計算               */
void moveParticle( void );                /* 速度・位置の仮更新         */
void collision( void );                   /* 粒子衝突処理               */
void calPressure( void );                 /* 圧力計算 (全体の司令塔)    */
void calNumberDensity( void );            /* 粒子数密度 n の計算        */
void setBoundaryCondition( void );        /* 自由表面/内部の判定        */
void setSourceTerm( void );               /* 連立一次方程式の右辺 {b}   */
void setMatrix( void );                   /* 連立一次方程式の行列 [A]   */
void exceptionalProcessingForBoundaryCondition( void ); /* 境界条件の例外処理 */
void checkBoundaryCondition( void );      /* ディリクレ境界の連結性検査 */
void increaseDiagonalTerm( void );        /* 対角項を 2 倍にする救済処置 */
void solveSimultaniousEquationsByGaussEliminationMethod( void );
                                          /* ガウスの消去法で求解       */
void removeNegativePressure( void );      /* 負圧を 0 にする            */
void setMinimumPressure( void );          /* 近傍の最小圧力を記録       */
void calPressureGradient( void );         /* 圧力勾配 → 加速度          */
void moveParticleUsingPressureGradient( void ); /* 速度・位置の修正     */
void writeData_inProfFormat( void );      /* prof 形式で出力            */
void writeData_inVtuFormat( void );       /* VTK(vtu) 形式で出力        */

/*=====================================================================
  【グローバル変数 (粒子ごとの配列)】

  ●配列の添字の決まり
    ・スカラー量 (圧力など)   … 粒子番号 i を使って  X[i]
    ・ベクトル量 (位置など)   … X[i*3  ] が x 成分
                                X[i*3+1] が y 成分
                                X[i*3+2] が z 成分
      2 次元計算でも z 成分の領域は確保しており、値は常に 0 である。

  ●static を付けているのは、このファイル内だけで使う変数であることを
    明示するため (他のファイルから参照できないようにする)。
=====================================================================*/
static double Acceleration[3*ARRAY_SIZE];
                        /* 各粒子の加速度 [m/s^2]。
                           重力・粘性・圧力勾配をここに足し込み、
                           速度更新に使ったあと 0 にクリアされる     */
static int    ParticleType[ARRAY_SIZE];
                        /* 各粒子の種類。FLUID / WALL / DUMMY_WALL /
                           GHOST のいずれかが入る                   */
static double Position[3*ARRAY_SIZE];
                        /* 各粒子の位置ベクトル [m]                 */
static double Velocity[3*ARRAY_SIZE];
                        /* 各粒子の速度ベクトル [m/s]               */
static double Pressure[ARRAY_SIZE];
                        /* 各粒子の圧力 [Pa]。
                           圧力ポアソン方程式の解 (未知数)          */
static double NumberDensity[ARRAY_SIZE];
                        /* 各粒子の粒子数密度 n [-]。
                           近傍粒子の重み w の総和。流体の密度に比例
                           する量で、自由表面判定にも使う           */
static int    BoundaryCondition[ARRAY_SIZE];
                        /* 圧力計算上の分類。
                           INNER_PARTICLE / SURFACE_PARTICLE /
                           GHOST_OR_DUMMY のいずれか               */
static double SourceTerm[ARRAY_SIZE];
                        /* 圧力ポアソン方程式の右辺ベクトル {b}。
                           ガウスの消去法の途中で書き換えられる     */
static int    FlagForCheckingBoundaryCondition[ARRAY_SIZE];
                        /* ディリクレ境界の連結性チェック用の作業配列 */
static double CoefficientMatrix[ARRAY_SIZE * ARRAY_SIZE];
                        /* 圧力ポアソン方程式の係数行列 [A]。
                           1 次元配列に 2 次元行列を詰めており、
                           i 行 j 列は CoefficientMatrix[i*n+j]
                           (n = NumberOfParticles) でアクセスする。
                           約 200MB を消費する巨大配列              */
static double MinimumPressure[ARRAY_SIZE];
                        /* 各粒子の近傍における最小圧力 [Pa]。
                           圧力勾配モデルを安定化させる (粒子間に
                           必ず斥力が働くようにする) ために使う     */

int    FileNumber;      /* 出力ファイルの通し番号。出力のたびに +1   */
double Time;            /* 現在の計算時刻 [s]                       */
int    NumberOfParticles;
                        /* 全粒子数。初期配置関数で決まる           */

double Re_forNumberDensity,  Re2_forNumberDensity;
                        /* 粒子数密度用の影響半径 re とその 2 乗    */
double Re_forGradient,       Re2_forGradient;
                        /* 勾配モデル用の影響半径 re とその 2 乗    */
double Re_forLaplacian,      Re2_forLaplacian;
                        /* ラプラシアンモデル用の影響半径とその2乗  */
                        /* 2 乗を持っておくのは、距離比較のたびに
                           sqrt() を呼ばずに済ませるための高速化    */

double N0_forNumberDensity;
                        /* 基準粒子数密度 n0 (粒子数密度用 re)。
                           粒子が初期格子状に隙間なく詰まった理想
                           状態での n の値。非圧縮条件 n = n0 の基準 */
double N0_forGradient;  /* 基準粒子数密度 n0 (勾配用 re)            */
double N0_forLaplacian; /* 基準粒子数密度 n0 (ラプラシアン用 re)    */
double Lambda;          /* ラプラシアンモデルの係数 λ。
                           λ = Σ(r^2・w) / Σw で定義され、
                           拡散方程式の解と統計的に一致させるための
                           補正係数                                 */
double collisionDistance, collisionDistance2;
                        /* 衝突判定距離とその 2 乗                  */
double FluidDensity;    /* 流体密度 ρ [kg/m^3]                      */


/*=====================================================================
  【関数名】main
  【機能】  プログラムの入口。粒子の初期配置 → 定数計算 →
            時間発展ループ、という全体の流れを制御する。
  【引数】  int    argc : コマンドライン引数の個数 (未使用)
            char** argv : コマンドライン引数の文字列配列 (未使用)
  【戻り値】int : 正常終了なら 0 を返す (OS に渡される終了コード)
  【呼び出し元】C ランタイム (プログラム起動時に自動的に呼ばれる)
=====================================================================*/
int main( int argc, char** argv ) {

  /* 本プログラムはコマンドライン引数を使わない。
     VC++ の警告レベル /W4 で出る「未使用の引数」警告 C4100 を
     避けるため、明示的に void へキャストして「使わない」ことを示す */
  (void)argc;
  (void)argv;

  printf("\n*** START PARTICLE-SIMULATION ***\n");

  /* 次元数 DIM に応じて初期配置関数を切り替える。
     DIM はマクロ定数なので、実際にはコンパイル時にどちらか一方が
     選ばれる                                                        */
  if( DIM == 2 ){
    initializeParticlePositionAndVelocity_for2dim();
  }else{
    initializeParticlePositionAndVelocity_for3dim();
  }

  calConstantParameter();   /* 影響半径・n0・λ など計算中不変の量を求める */
  mainLoopOfSimulation();   /* ここで時間発展の全ステップを回す          */

  printf("*** END ***\n\n");
  return 0;
}


/*=====================================================================
  【関数名】initializeParticlePositionAndVelocity_for2dim
  【機能】  2 次元ダムブレイク問題の初期粒子配置を作る。
            格子状に候補点を走査し、その点がどの領域 (ダミー壁／壁／
            空／流体) に属するかを判定して粒子を生成する。
            生成した粒子の種類・位置を設定し、速度を 0 に初期化する。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用 (書き換えるグローバル変数)】
            ParticleType[]    : 各粒子の種類
            Position[]        : 各粒子の初期位置
            Velocity[]        : 全成分 0 に初期化
            NumberOfParticles : 生成された粒子の総数
  【呼び出し元】main()  (DIM == 2 のとき)
  【計算領域】
            水槽 : x = 0.0〜1.0 [m], y = 0.0〜0.6 [m]
            水柱 : x = 0.0〜0.25[m], y = 0.0〜0.50[m]
            壁は水槽の外側に 2 層、さらにその外にダミー壁を配置する。
            領域判定は「後に書いた if 文が優先される」ように、
            ダミー壁 → 壁 → 空 → 流体 の順に上書きしている。
=====================================================================*/
void initializeParticlePositionAndVelocity_for2dim( void ){
  int iX, iY;      /* 候補点の格子番号 (x 方向, y 方向)              */
  int nX, nY;      /* 格子番号の上限値 (走査する格子点の個数)        */
  double x, y, z;  /* 候補点の座標 [m]。2 次元なので z は常に 0.0    */
  int i = 0;       /* 生成済み粒子の個数 = 次に格納する粒子番号      */
  int flagOfParticleGeneration; /* この候補点に粒子を作るか
                                   (ON:作る / OFF:作らない)          */

  /* 走査する格子点の範囲を決める。
     +5 と、下の -4 から始める分で、水槽の外側 4 層まで走査できる    */
  nX = (int)(1.0/PARTICLE_DISTANCE)+5;
  nY = (int)(0.6/PARTICLE_DISTANCE)+5;

  for(iX= -4;iX<nX;iX++){
    for(iY= -4;iY<nY;iY++){
      /* 格子番号から実座標へ変換する */
      x = PARTICLE_DISTANCE * (double)(iX);
      y = PARTICLE_DISTANCE * (double)(iY);
      z = 0.0;
      flagOfParticleGeneration = OFF;  /* まずは「作らない」で初期化 */

      /* dummy wall region */
      /* ダミー壁領域 : 水槽の外側 4 層分。
         壁粒子の粒子数密度が不足しないように埋めておくための粒子   */
      if( ((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-4.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)) ){
	ParticleType[i]=DUMMY_WALL;
	flagOfParticleGeneration = ON;
      }

      /* wall region */
      /* 壁領域 : 水槽の外側 2 層分。ダミー壁を上書きする            */
      if( ((x>-2.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+2.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)) ){
	ParticleType[i]=WALL;
	flagOfParticleGeneration = ON;
      }

      /* wall region */
      /* 壁領域 (上部) : y = 0.6 付近の天井部分を壁にする            */
      if( ((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.6-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)) ){
	ParticleType[i]=WALL;
	flagOfParticleGeneration = ON;
      }

      /* empty region */
      /* 空領域 : 水槽の内側は空気なので粒子を作らない。
         上で壁と判定された内部の点をここで取り消している           */
      if( ((x>0.0+EPS)&&(x<=1.00+EPS))&&( y>0.0+EPS )){
	flagOfParticleGeneration = OFF;
      }

      /* fluid region */
      /* 流体領域 : 左側の水柱 (幅 0.25m × 高さ 0.50m)              */
      if( ((x>0.0+EPS)&&(x<=0.25+EPS)) &&((y>0.0+EPS)&&(y<=0.50+EPS)) ){
	ParticleType[i]=FLUID;
	flagOfParticleGeneration = ON;
      }

      /* 粒子を作ると決まったら、位置を格納して粒子番号を 1 進める  */
      if( flagOfParticleGeneration == ON){
	Position[i*3]=x; Position[i*3+1]=y; Position[i*3+2]=z;
	i++;
      }
    }
  }

  NumberOfParticles = i;   /* 走査終了時点の i が総粒子数になる      */

  /* 全粒子の速度を 0 に初期化する (静止状態からスタート)。
     ベクトル量なので粒子数の 3 倍の要素をまとめてクリアする        */
  for(i=0;i<NumberOfParticles*3;i++) { Velocity[i]=0.0; }
}


/*=====================================================================
  【関数名】initializeParticlePositionAndVelocity_for3dim
  【機能】  3 次元ダムブレイク問題の初期粒子配置を作る。
            2 次元版に z 方向 (奥行き 0.3m) の判定を加えたもので、
            処理の考え方はまったく同じである。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用 (書き換えるグローバル変数)】
            ParticleType[], Position[], Velocity[], NumberOfParticles
  【呼び出し元】main()  (DIM != 2 のとき)
  【注意】  3 次元では粒子数が一気に増えるため、ARRAY_SIZE と
            PARTICLE_DISTANCE の設定に注意すること。
=====================================================================*/
void initializeParticlePositionAndVelocity_for3dim( void ){
  int iX, iY, iZ;   /* 候補点の格子番号 (x, y, z 方向)               */
  int nX, nY, nZ;   /* 各方向の格子番号の上限値                      */
  double x, y, z;   /* 候補点の座標 [m]                              */
  int i = 0;        /* 生成済み粒子の個数 = 次に格納する粒子番号     */
  int flagOfParticleGeneration;  /* 粒子を作るかどうかのフラグ       */

  nX = (int)(1.0/PARTICLE_DISTANCE)+5;
  nY = (int)(0.6/PARTICLE_DISTANCE)+5;
  nZ = (int)(0.3/PARTICLE_DISTANCE)+5;

  for(iX= -4;iX<nX;iX++){
    for(iY= -4;iY<nY;iY++){
      for(iZ= -4;iZ<nZ;iZ++){
	x = PARTICLE_DISTANCE * iX;
	y = PARTICLE_DISTANCE * iY;
	z = PARTICLE_DISTANCE * iZ;
	flagOfParticleGeneration = OFF;

	/* dummy wall region */
	/* ダミー壁領域 : 計算領域の外側 4 層分 (x,y,z すべて)        */
	if( (((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-4.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-4.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+4.0*PARTICLE_DISTANCE+EPS ))){
	  ParticleType[i]=DUMMY_WALL;
	  flagOfParticleGeneration = ON;
	}

	/* wall region */
	/* 壁領域 : 計算領域の外側 2 層分                             */
	if( (((x>-2.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+2.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-2.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+2.0*PARTICLE_DISTANCE+EPS ))){
	  ParticleType[i]=WALL;
	  flagOfParticleGeneration = ON;
	}

	/* wall region */
	/* 壁領域 (上部) : y = 0.6 付近の天井部分                     */
	if( (((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.6-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-4.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+4.0*PARTICLE_DISTANCE+EPS ))){
	  ParticleType[i]=WALL;
	  flagOfParticleGeneration = ON;
	}

	/* empty region */
	/* 空領域 : 水槽の内側 (空気) には粒子を作らない              */
	if( (((x>0.0+EPS)&&(x<=1.00+EPS))&&( y>0.0+EPS ))&&( (z>0.0+EPS )&&(z<=0.3+EPS ))){
	  flagOfParticleGeneration = OFF;
	}

	/* fluid region */
	/* 流体領域 : 水柱 (0.25m × 0.5m × 0.3m)                     */
	if( (((x>0.0+EPS)&&(x<=0.25+EPS))&&( (y>0.0+EPS)&&(y<0.5+EPS) ))&&( (z>0.0+EPS )&&(z<=0.3+EPS ))){
	  ParticleType[i]=FLUID;
	  flagOfParticleGeneration = ON;
	}

	if( flagOfParticleGeneration == ON){
	  Position[i*3  ]=x;
	  Position[i*3+1]=y;
	  Position[i*3+2]=z;
	  i++;
	}
      }
    }
  }

  NumberOfParticles = i;   /* 総粒子数を確定                         */

  /* 全粒子の速度を 0 に初期化する                                  */
  for(i=0;i<NumberOfParticles*3;i++) { Velocity[i]=0.0; }
}


/*=====================================================================
  【関数名】calConstantParameter
  【機能】  計算中に変化しない定数を、マクロ定数からグローバル変数へ
            代入し、必要な派生量 (2 乗値、n0、λ) を求めておく。
            毎ステップ再計算しないで済ませるための前処理である。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用 (書き換えるグローバル変数)】
            Re_forNumberDensity / Re_forGradient / Re_forLaplacian
            Re2_* (それぞれの 2 乗)
            FluidDensity, collisionDistance, collisionDistance2
            FileNumber, Time
            さらに calNZeroAndLambda() 経由で N0_* と Lambda
  【呼び出し元】main()  (初期配置の直後に 1 回だけ)
=====================================================================*/
void calConstantParameter( void ){

  /* 3 種類の影響半径をマクロ定数から取り込む                       */
  Re_forNumberDensity  = RADIUS_FOR_NUMBER_DENSITY;
  Re_forGradient       = RADIUS_FOR_GRADIENT;
  Re_forLaplacian      = RADIUS_FOR_LAPLACIAN;

  /* 距離の比較では sqrt() を避けて 2 乗同士で比較したいので、
     影響半径の 2 乗もあらかじめ計算しておく                        */
  Re2_forNumberDensity = Re_forNumberDensity*Re_forNumberDensity;
  Re2_forGradient      = Re_forGradient*Re_forGradient;
  Re2_forLaplacian     = Re_forLaplacian*Re_forLaplacian;

  calNZeroAndLambda();  /* 基準粒子数密度 n0 とラプラシアン係数 λ   */

  FluidDensity       = FLUID_DENSITY;
  collisionDistance  = COLLISION_DISTANCE;
  collisionDistance2 = collisionDistance*collisionDistance;

  FileNumber=0;   /* 出力ファイル番号を 0 から始める                */
  Time=0.0;       /* 計算開始時刻                                    */
}


/*=====================================================================
  【関数名】calNZeroAndLambda
  【機能】  基準粒子数密度 n0 とラプラシアンモデルの係数 λ を求める。

            n0 : 粒子が初期格子状に無限に詰まっている理想状態で、
                 ある 1 粒子が感じる粒子数密度。
                 原点に仮想の粒子 i を置き、その周囲の格子点 j に
                 粒子があるものとして重み w の総和をとる。
                 MPS 法では「n = n0 を保つこと」が非圧縮条件になる。

            λ  : ラプラシアンモデルの補正係数。
                 λ = Σ(|rij|^2・w) / Σw
                 粒子法のラプラシアン近似が拡散方程式の解と
                 統計的に一致するようにするための係数。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用 (書き換えるグローバル変数)】
            N0_forNumberDensity, N0_forGradient, N0_forLaplacian, Lambda
  【呼び出し元】calConstantParameter()
=====================================================================*/
void calNZeroAndLambda( void ){
  int iX, iY, iZ;          /* 仮想格子点の番号                       */
  int iZ_start, iZ_end;    /* z 方向の走査範囲 (次元によって変わる)  */
  double xj, yj, zj;       /* 近傍粒子 j (仮想格子点) の座標         */
  double distance, distance2; /* 粒子 i-j 間の距離とその 2 乗        */
  double xi, yi, zi;       /* 着目粒子 i の座標 (原点に置く)         */

  /* 2 次元なら z 方向は 1 層だけ (iZ = 0) 走査する。
     3 次元なら -4〜4 の 9 層を走査する                             */
  if( DIM == 2 ){
    iZ_start = 0; iZ_end = 1;
  }else{
    iZ_start = -4; iZ_end = 5;
  }

  /* 総和をとる変数をゼロクリアしてから積算する                     */
  N0_forNumberDensity = 0.0;
  N0_forGradient      = 0.0;
  N0_forLaplacian     = 0.0;
  Lambda              = 0.0;

  xi = 0.0;  yi = 0.0;  zi = 0.0;   /* 着目粒子は原点に置く          */

  for(iX= -4;iX<5;iX++){
    for(iY= -4;iY<5;iY++){
      for(iZ= iZ_start;iZ<iZ_end;iZ++){
	/* 自分自身 (原点) は和に含めないので飛ばす                  */
	if( ((iX==0)&&(iY==0)) && (iZ==0) )continue;

	/* 格子点の座標を求める                                      */
	xj = PARTICLE_DISTANCE * (double)(iX);
	yj = PARTICLE_DISTANCE * (double)(iY);
	zj = PARTICLE_DISTANCE * (double)(iZ);

	/* 粒子間距離を計算する                                      */
	distance2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi);
	distance = sqrt(distance2);

	/* 影響半径ごとに重みを積算する。
	   distance >= re の格子点は weight() が 0 を返すので、
	   自動的に和から除外される                                  */
	N0_forNumberDensity += weight(distance, Re_forNumberDensity);
	N0_forGradient      += weight(distance, Re_forGradient);
	N0_forLaplacian     += weight(distance, Re_forLaplacian);

	/* λ の分子 Σ(r^2・w) を積算する                            */
	Lambda              += distance2 * weight(distance, Re_forLaplacian);
      }
    }
  }

  /* λ = Σ(r^2・w) / Σw  (分母は N0_forLaplacian そのもの)         */
  Lambda = Lambda/N0_forLaplacian;
}


/*=====================================================================
  【関数名】weight
  【機能】  MPS 法の重み関数 (カーネル関数) w(r, re) を計算する。

                     ┌ re/r − 1   (r < re)
              w(r) = │
                     └ 0          (r >= re)

            距離が近いほど大きな値をとり、影響半径 re で 0 になる。
            r → 0 で発散するため、粒子どうしが極端に接近しても
            斥力が働き、粒子の重なりを防ぐ効果がある。
  【引数】  double distance : 粒子間の距離 r [m]
            double re       : 影響半径 re [m]
  【戻り値】double : 重み w の値 [-] (0 以上)
  【呼び出し元】calNZeroAndLambda(), calViscosity(), calNumberDensity(),
                setMatrix(), calPressureGradient()
  【注意】  distance == 0.0 でこの関数を呼ぶとゼロ除算になる。
            呼び出し側で必ず j != i としているため、実際には
            distance が 0 になることはない。
=====================================================================*/
double weight( double distance, double re ){
  double weightIJ;   /* 粒子 i-j 間の重みの値 (戻り値になる)        */

  if( distance >= re ){
    weightIJ = 0.0;                   /* 影響半径の外 → 相互作用なし */
  }else{
    weightIJ = (re/distance) - 1.0;   /* 影響半径の内側              */
  }
  return weightIJ;
}


/*=====================================================================
  【関数名】mainLoopOfSimulation
  【機能】  時間発展の main ループ。FINISH_TIME に達するまで、
            1 ステップ分の計算 (陽的計算 → 圧力計算 → 修正) を
            繰り返し実行し、一定間隔で結果をファイルに出力する。

            ★この関数を読めば MPS 法 1 ステップの手順が分かる★

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Time を DT ずつ進め、粒子の状態量すべてを更新する。
            また結果ファイルを書き出す。
  【呼び出し元】main()
=====================================================================*/
void mainLoopOfSimulation( void ){
  int iTimeStep = 0;   /* 実行済みの時間ステップ数 (カウンタ)       */

  /* 計算開始前に初期状態 (t = 0) を出力しておく                    */
  writeData_inVtuFormat();
  writeData_inProfFormat();

  while(1){   /* 終了条件は while の内側で判定して break する        */

    /*--- (1) 陽的計算 : 圧力以外の力で仮の速度・位置を求める ---*/
    calGravity();    /* 重力による加速度を代入 (加速度を初期化も兼ねる) */
    calViscosity();  /* 粘性による加速度を加算                        */
    moveParticle();  /* 仮の速度 u* と仮の位置 r* に更新              */
    collision();     /* 接近しすぎた粒子どうしを反発させる            */

    /*--- (2) 陰的計算 : 非圧縮性を満たす圧力を求める ---*/
    calPressure();   /* 圧力ポアソン方程式を解いて Pressure[] を決定  */

    /*--- (3) 修正 : 圧力勾配で速度・位置を補正して確定させる ---*/
    calPressureGradient();              /* 圧力勾配 → 加速度         */
    moveParticleUsingPressureGradient();/* 速度 u^(k+1)・位置 r^(k+1) */

    /*--- 時刻を進める ---*/
    iTimeStep++;
    Time += DT;

    /* OUTPUT_INTERVAL ステップごとに経過表示とファイル出力を行う    */
    if( (iTimeStep % OUTPUT_INTERVAL) == 0 ){
      printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticless: %d\n", iTimeStep, Time, NumberOfParticles);
      writeData_inVtuFormat();
      writeData_inProfFormat();
    }

    /* 終了時刻に達したらループを抜ける                              */
    if( Time >= FINISH_TIME ){break;}
  }
}


/*=====================================================================
  【関数名】calGravity
  【機能】  全粒子の加速度に重力加速度を「代入」する。
            加算ではなく代入なので、この関数はステップ先頭で
            加速度配列を初期化する役割も兼ねている。
            流体粒子だけが重力を受け、壁粒子などは 0 のままとする
            (壁は動かさないため)。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Acceleration[] を書き換える
  【呼び出し元】mainLoopOfSimulation()  (1 ステップの最初)
=====================================================================*/
void calGravity( void ){
  int i;   /* 粒子番号のループカウンタ                              */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      /* 流体粒子には重力加速度 (G_X, G_Y, G_Z) を与える            */
      Acceleration[i*3  ]=G_X;
      Acceleration[i*3+1]=G_Y;
      Acceleration[i*3+2]=G_Z;
    }else{
      /* 流体以外 (壁・ダミー壁・ゴースト) は加速度 0               */
      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;
    }
  }
}


/*=====================================================================
  【関数名】calViscosity
  【機能】  粘性項 ν∇^2 u を MPS 法のラプラシアンモデルで離散化し、
            加速度に加算する。

            ラプラシアンモデル :
              <∇^2 φ>i = (2・DIM)/(λ・n0) ・ Σ (φj − φi)・w(rij)

            これを速度 u に適用し、動粘性係数 ν を掛けたものが
            粘性による加速度になる。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Acceleration[] に粘性項を加算する
  【呼び出し元】mainLoopOfSimulation()  (calGravity の直後)
  【計算量】全粒子ペアを総当たりするので O(N^2)
=====================================================================*/
void calViscosity( void ){
  int i,j;   /* i:着目粒子の番号, j:近傍粒子の番号                   */
  double viscosityTerm_x, viscosityTerm_y, viscosityTerm_z;
             /* 粒子 i が受ける粘性項 (近傍の総和を貯める変数)       */
  double distance, distance2;  /* 粒子 i-j 間の距離とその 2 乗       */
  double w;                    /* 粒子 i-j 間の重み w(rij)           */
  double xij, yij, zij;        /* 相対位置ベクトル rij = rj − ri     */
  double a;                    /* ラプラシアンモデルの係数
                                  a = ν(2・DIM)/(n0・λ)             */

  /* 係数は粒子によらず一定なので、ループの外で 1 回だけ計算する    */
  a = (KINEMATIC_VISCOSITY)*(2.0*DIM)/(N0_forLaplacian*Lambda);

  for(i=0;i<NumberOfParticles;i++){
    /* 粘性で動かすのは流体粒子だけなので、それ以外は飛ばす         */
    if(ParticleType[i] != FLUID) continue;

    /* 総和用の変数をゼロクリアする                                 */
    viscosityTerm_x = 0.0;  viscosityTerm_y = 0.0;  viscosityTerm_z = 0.0;

    for(j=0;j<NumberOfParticles;j++){
      /* 自分自身と計算対象外の粒子は除外する                       */
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;

      /* 相対位置ベクトルと距離を求める                             */
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);

      /* 影響半径の内側にある粒子だけを近傍粒子として扱う           */
      if(distance<Re_forLaplacian){
	w =  weight(distance, Re_forLaplacian);
	/* Σ(uj − ui)・w を成分ごとに積算する                       */
	viscosityTerm_x +=(Velocity[j*3  ]-Velocity[i*3  ])*w;
	viscosityTerm_y +=(Velocity[j*3+1]-Velocity[i*3+1])*w;
	viscosityTerm_z +=(Velocity[j*3+2]-Velocity[i*3+2])*w;
      }
    }

    /* 総和に係数 a を掛けて粘性による加速度にする                  */
    viscosityTerm_x = viscosityTerm_x * a;
    viscosityTerm_y = viscosityTerm_y * a;
    viscosityTerm_z = viscosityTerm_z * a;

    /* 重力の加速度に加算する (代入ではなく += であることに注意)    */
    Acceleration[i*3  ] += viscosityTerm_x;
    Acceleration[i*3+1] += viscosityTerm_y;
    Acceleration[i*3+2] += viscosityTerm_z;
  }
}


/*=====================================================================
  【関数名】moveParticle
  【機能】  重力・粘性で求めた加速度を使って、流体粒子の速度と位置を
            オイラー陽解法で更新する。ここで得られるのは圧力を考慮
            していない「仮の速度 u*」「仮の位置 r*」である。

              u* = u^k + a・Δt
              r* = r^k + u*・Δt

            更新後は次の計算のために加速度を 0 にクリアする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Velocity[], Position[] を更新し、Acceleration[] を 0 にする
  【呼び出し元】mainLoopOfSimulation()  (calViscosity の直後)
=====================================================================*/
void moveParticle( void ){
  int i;   /* 粒子番号のループカウンタ                              */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      /* 速度を更新 : u* = u + a・Δt                                */
      Velocity[i*3  ] += Acceleration[i*3  ]*DT;
      Velocity[i*3+1] += Acceleration[i*3+1]*DT;
      Velocity[i*3+2] += Acceleration[i*3+2]*DT;

      /* 位置を更新 : r* = r + u*・Δt  (更新後の速度を使う)         */
      Position[i*3  ] += Velocity[i*3  ]*DT;
      Position[i*3+1] += Velocity[i*3+1]*DT;
      Position[i*3+2] += Velocity[i*3+2]*DT;
    }
    /* 流体か否かによらず、使い終わった加速度をクリアしておく       */
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
  }
}


/*=====================================================================
  【関数名】collision
  【機能】  粒子どうしが異常に接近 (めり込み) した場合に、剛体球の
            衝突として運動量保存則に基づき速度を補正する。
            圧力計算だけでは防ぎきれない粒子の重なりを防止し、
            計算の破綻を避けるための安定化処理である。

            距離が collisionDistance より小さく、かつ互いに近づく
            向きに動いている粒子ペアに対して、力積

              forceDT = (1+e)・mi・mj/(mi+mj) ・ (相対速度の法線成分)

            を与え、速度を跳ね返す。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Velocity[], Position[] を補正する
  【呼び出し元】mainLoopOfSimulation()  (moveParticle の直後)
  【重要】  速度の補正は一時配列 VelocityAfterCollision[] に貯めてから
            最後にまとめて反映する。計算途中で Velocity[] を書き換えて
            しまうと、粒子を処理する順番によって結果が変わってしまう
            ためである。
=====================================================================*/
void collision( void ){
  int    i,j;              /* i:着目粒子, j:相手粒子                 */
  double xij, yij, zij;    /* 相対位置ベクトル rij = rj − ri         */
  double distance,distance2; /* 粒子間距離とその 2 乗                */
  double forceDT; /* forceDT is the impulse of collision between particles */
                  /* 衝突によって生じる力積 (力 × Δt)。
                     最初は相対速度の法線成分を入れ、後で係数を掛けて
                     力積に変換している                              */
  double mi, mj;  /* 粒子 i, j の質量に相当する量。
                     本プログラムでは体積が等しいので密度で代用する  */
  double velocity_ix, velocity_iy, velocity_iz;
                  /* 粒子 i の速度の作業用コピー。
                     j のループ内で順次更新していく                  */
  double e = COEFFICIENT_OF_RESTITUTION;   /* 反発係数               */
  static double VelocityAfterCollision[3*ARRAY_SIZE];
                  /* 衝突処理後の速度を一時的に保持する配列。
                     static にしているのはスタックオーバーフローを
                     避けるため (3×5000×8byte = 120KB)              */

  /* まず現在の速度をそのままコピーしておく                         */
  for(i=0;i<3*NumberOfParticles;i++){
    VelocityAfterCollision[i] = Velocity[i];
  }

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      mi = FluidDensity;                /* 粒子 i の質量 (密度で代用) */
      velocity_ix = Velocity[i*3  ];    /* 粒子 i の速度を作業変数へ  */
      velocity_iy = Velocity[i*3+1];
      velocity_iz = Velocity[i*3+2];

      for(j=0;j<NumberOfParticles;j++){
	if( (j==i) || (ParticleType[j]==GHOST) ) continue;

	/* 相対位置と距離の 2 乗を求める (sqrt はまだ計算しない)     */
	xij = Position[j*3  ] - Position[i*3  ];
	yij = Position[j*3+1] - Position[i*3+1];
	zij = Position[j*3+2] - Position[i*3+2];
	distance2 = (xij*xij) + (yij*yij) + (zij*zij);

	/* 衝突判定距離より近い場合だけ衝突処理を行う                */
	if(distance2<collisionDistance2){
	  distance = sqrt(distance2);

	  /* 相対速度 (ui − uj) の、粒子 i から j に向かう単位ベクトル
	     方向の成分 (法線成分) を求める。
	     正なら「近づいている」= 衝突している                    */
	  forceDT = (velocity_ix-Velocity[j*3  ])*(xij/distance)
	           +(velocity_iy-Velocity[j*3+1])*(yij/distance)
	           +(velocity_iz-Velocity[j*3+2])*(zij/distance);

	  if(forceDT > 0.0){   /* 近づいているときだけ跳ね返す        */
	    mj = FluidDensity; /* 相手の質量 (壁も同じ値として扱う)   */

	    /* 換算質量と反発係数から力積を求める                     */
	    forceDT *= (1.0+e)*mi*mj/(mi+mj);

	    /* 力積 ÷ 質量 = 速度変化。j から遠ざかる向きに減速する   */
	    velocity_ix -= (forceDT/mi)*(xij/distance);
	    velocity_iy -= (forceDT/mi)*(yij/distance);
	    velocity_iz -= (forceDT/mi)*(zij/distance);

	    /* 衝突が起きたことを知りたい場合は下のコメントを外す。
	       ただし大量に出力されるので通常は無効にしておく         */
	    /*
	    if(j>i){ fprintf(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
	    */
	  }
	}
      }
      /* 粒子 i について求めた衝突後速度を一時配列へ格納する         */
      VelocityAfterCollision[i*3  ] = velocity_ix;
      VelocityAfterCollision[i*3+1] = velocity_iy;
      VelocityAfterCollision[i*3+2] = velocity_iz;
    }
  }

  /* すべての粒子の衝突計算が終わってから、まとめて反映する          */
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      /* 速度変化分だけ位置も補正する (Δu・Δt だけずらす)           */
      Position[i*3  ] += (VelocityAfterCollision[i*3  ]-Velocity[i*3  ])*DT;
      Position[i*3+1] += (VelocityAfterCollision[i*3+1]-Velocity[i*3+1])*DT;
      Position[i*3+2] += (VelocityAfterCollision[i*3+2]-Velocity[i*3+2])*DT;

      /* 速度を衝突後の値で置き換える                                */
      Velocity[i*3  ] = VelocityAfterCollision[i*3  ];
      Velocity[i*3+1] = VelocityAfterCollision[i*3+1];
      Velocity[i*3+2] = VelocityAfterCollision[i*3+2];
    }
  }
}


/*=====================================================================
  【関数名】calPressure
  【機能】  圧力を求める一連の処理をまとめて呼び出す「司令塔」関数。
            この関数自体は計算をせず、以下の 7 つの関数を決まった順に
            呼ぶだけである。

              1. calNumberDensity()   粒子数密度 n を数える
              2. setBoundaryCondition() 自由表面か内部かを判定
              3. setSourceTerm()      連立方程式の右辺 {b} を作る
              4. setMatrix()          連立方程式の係数行列 [A] を作る
              5. solveSimultanious…  [A]{P}={b} をガウス消去法で解く
              6. removeNegativePressure() 負の圧力を 0 にする
              7. setMinimumPressure() 近傍最小圧力を記録する

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Pressure[], MinimumPressure[] などを更新する
  【呼び出し元】mainLoopOfSimulation()  (collision の直後)
=====================================================================*/
void calPressure( void ){
  calNumberDensity();     /* 粒子数密度 n を計算する                 */
  setBoundaryCondition(); /* 自由表面粒子 / 内部粒子を判定する       */
  setSourceTerm();        /* 右辺ベクトル {b} を作る                 */
  setMatrix();            /* 係数行列 [A] を作る                     */
  solveSimultaniousEquationsByGaussEliminationMethod(); /* 解く       */
  removeNegativePressure(); /* 負圧をカットする                      */
  setMinimumPressure();   /* 勾配計算用に近傍最小圧力を求める        */
}


/*=====================================================================
  【関数名】calNumberDensity
  【機能】  各粒子の粒子数密度 n を計算する。

              ni = Σ_{j≠i} w(|rj − ri|, re)

            n は流体の密度に比例する量で、MPS 法では
              ・非圧縮条件 (n を n0 に保つ)
              ・自由表面の判定 (n が小さい粒子は表面)
            の両方に使われる、最も基本的な量である。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】NumberDensity[] を更新する
  【呼び出し元】calPressure()
  【計算量】O(N^2)
=====================================================================*/
void calNumberDensity( void ){
  int    i,j;             /* i:着目粒子, j:近傍粒子                  */
  double xij, yij, zij;   /* 相対位置ベクトル                        */
  double distance, distance2; /* 粒子間距離とその 2 乗               */
  double w;               /* 重み w(rij)                             */

  for(i=0;i<NumberOfParticles;i++){
    NumberDensity[i] = 0.0;   /* 総和の前にゼロクリア                */
    if(ParticleType[i] == GHOST) continue;  /* 計算対象外            */

    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);

      /* 影響半径の外なら weight() が 0 を返すので、
	 if 文で距離判定をしなくても結果は正しくなる               */
      w =  weight(distance, Re_forNumberDensity);
      NumberDensity[i] += w;
    }
  }
}


/*=====================================================================
  【関数名】setBoundaryCondition
  【機能】  各粒子を圧力計算上の 3 種類に分類する。

              GHOST_OR_DUMMY   : 圧力を解かない粒子
              SURFACE_PARTICLE : 自由表面粒子 → 圧力 P = 0 を与える
                                 (ディリクレ境界条件)
              INNER_PARTICLE   : 内部粒子 → 圧力を未知数として解く

            自由表面の判定は粒子数密度で行う。表面付近では片側に
            近傍粒子が存在しないため n が小さくなる性質を利用し、
              n < β・n0   (β = THRESHOLD_RATIO_OF_NUMBER_DENSITY)
            なら表面とみなす。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】BoundaryCondition[] を更新する
  【呼び出し元】calPressure()
=====================================================================*/
void setBoundaryCondition( void ){
  int i;                              /* 粒子番号                    */
  double n0 = N0_forNumberDensity;    /* 基準粒子数密度              */
  double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY; /* 判定しきい値 β */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ){
      /* ゴースト粒子とダミー壁粒子は圧力を解かない                 */
      BoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if( NumberDensity[i] < beta * n0 ){
      /* 粒子数密度が不足 → 自由表面 (P = 0 のディリクレ境界)       */
      BoundaryCondition[i]=SURFACE_PARTICLE;
    }else{
      /* 周囲が粒子で満たされている → 内部粒子 (未知数)             */
      BoundaryCondition[i]=INNER_PARTICLE;
    }
  }
}


/*=====================================================================
  【関数名】setSourceTerm
  【機能】  圧力ポアソン方程式 [A]{P} = {b} の右辺ベクトル {b} を作る。

            内部粒子では、粒子数密度 n を基準値 n0 に戻すための
            ソース項を与える :

              bi = γ ・ (1/Δt^2) ・ (ni − n0)/n0

            γ は緩和係数 (RELAXATION_COEFFICIENT_FOR_PRESSURE)。
            1 ステップで一気に補正すると圧力が振動するため、
            γ(=0.2) を掛けて緩やかに補正する。

            自由表面粒子は P = 0 を与えるので bi = 0 とする。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】SourceTerm[] を更新する
  【呼び出し元】calPressure()
=====================================================================*/
void setSourceTerm( void ){
  int i;                                        /* 粒子番号          */
  double n0    = N0_forNumberDensity;           /* 基準粒子数密度    */
  double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE; /* 緩和係数 γ  */

  for(i=0;i<NumberOfParticles;i++){
    SourceTerm[i]=0.0;   /* まず 0 で初期化する                      */

    /* 圧力を解かない粒子は 0 のままにしておく                       */
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ) continue;

    if(BoundaryCondition[i]==INNER_PARTICLE){
      /* 内部粒子 : 粒子数密度の n0 からのずれをソース項にする       */
      SourceTerm[i] = gamma * (1.0/(DT*DT))*((NumberDensity[i]-n0)/n0);
    }else if(BoundaryCondition[i]==SURFACE_PARTICLE){
      /* 自由表面粒子 : P = 0 (大気圧を基準) なので右辺も 0          */
      SourceTerm[i]=0.0;
    }
  }
}


/*=====================================================================
  【関数名】setMatrix
  【機能】  圧力ポアソン方程式 [A]{P} = {b} の係数行列 [A] を作る。

            ラプラシアンモデルを圧力に適用すると、粒子 i の行は

              対角項     A[i][i] = Σ_j (a・w(rij)/ρ) + 圧縮率項
              非対角項   A[i][j] = −a・w(rij)/ρ
              ただし     a = 2・DIM/(n0・λ)

            となる。対角項に COMPRESSIBILITY/Δt^2 を足すのは、
            わずかな圧縮性を許して行列を対角優位にし、計算を
            安定化させるためである。

            自由表面粒子の行は 0 のまま (P = 0 が既知) とし、
            最後に exceptionalProcessingForBoundaryCondition() を
            呼んで特異行列になっていないかを確認する。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】CoefficientMatrix[] を更新する
  【呼び出し元】calPressure()
  【計算量】行列のクリアと構築でそれぞれ O(N^2)
=====================================================================*/
void setMatrix( void ){
  double xij, yij, zij;       /* 相対位置ベクトル                    */
  double distance, distance2; /* 粒子間距離とその 2 乗               */
  double coefficientIJ;       /* 粒子 i-j ペアの係数                 */
  double n0 = N0_forLaplacian;/* ラプラシアン用の基準粒子数密度      */
  int    i,j;                 /* i:行番号, j:列番号 (どちらも粒子番号)*/
  double a;                   /* ラプラシアンモデルの係数
                                 a = 2・DIM/(n0・λ)                 */
  int n = NumberOfParticles;  /* 行列の一辺の大きさ。
                                 1 次元配列の添字計算 i*n+j に使う   */

  /* まず行列全体を 0 でクリアする                                  */
  for(i=0;i<NumberOfParticles;i++){
    for(j=0;j<NumberOfParticles;j++){
      CoefficientMatrix[i*n+j] = 0.0;
    }
  }

  a = 2.0*DIM/(n0*Lambda);

  for(i=0;i<NumberOfParticles;i++){
    /* 未知数として解くのは内部粒子だけ。
       それ以外の行は 0 のまま (後で対角項だけ扱う)                 */
    if(BoundaryCondition[i] != INNER_PARTICLE) continue;

    for(j=0;j<NumberOfParticles;j++){
      /* 自分自身と、圧力を解かない粒子は列に入れない               */
      if( (j==i) || (BoundaryCondition[j]==GHOST_OR_DUMMY) ) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);
      distance  = sqrt(distance2);

      /* 影響半径の外の粒子とは相互作用しない                       */
      if(distance>=Re_forLaplacian)continue;

      /* 粒子ペアごとの係数。密度 ρ で割って圧力→加速度の次元にする */
      coefficientIJ = a * weight(distance, Re_forLaplacian)/FluidDensity;

      CoefficientMatrix[i*n+j]  = (-1.0)*coefficientIJ; /* 非対角項  */
      CoefficientMatrix[i*n+i] += coefficientIJ;        /* 対角項に加算 */
    }
    /* 対角項に圧縮率の項を足して対角優位にし、解を安定化させる      */
    CoefficientMatrix[i*n+i] += (COMPRESSIBILITY)/(DT*DT);
  }

  /* ディリクレ境界が存在しない孤立粒子群への対策を行う              */
  exceptionalProcessingForBoundaryCondition();
}


/*=====================================================================
  【関数名】exceptionalProcessingForBoundaryCondition
  【機能】  ディリクレ境界条件 (自由表面の P = 0) を持たない粒子群が
            存在する場合の例外処理を行う。

            原文コメントの意図 :
              「流体にディリクレ境界条件が無い場合、例外処理として
                行列の対角項を増加させる。これによりディリクレ境界
                条件無しでも行列を解くことができる」

            圧力ポアソン方程式は、粒子のかたまりの中に P = 0 と
            定める点が 1 つも無いと解が一意に定まらない (行列が
            特異になり、ゼロ除算や発散を引き起こす)。
            そこで連結性を調べ、孤立している粒子の対角項を 2 倍に
            して無理やり解けるようにする。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】FlagForCheckingBoundaryCondition[] と
            CoefficientMatrix[] の対角項を書き換える
  【呼び出し元】setMatrix()  (行列を作り終えた直後)
=====================================================================*/
void exceptionalProcessingForBoundaryCondition( void ){
  /* If tere is no Dirichlet boundary condition on the fluid,
     increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
  checkBoundaryCondition();  /* 自由表面とつながっているかを調べる   */
  increaseDiagonalTerm();    /* つながっていない粒子の対角項を 2 倍  */
}


/*=====================================================================
  【関数名】checkBoundaryCondition
  【機能】  各粒子が「自由表面粒子 (ディリクレ境界) とつながって
            いるか」を幅優先探索の要領で調べる。

            手順 :
              1. 初期状態を設定する
                   圧力を解かない粒子 → GHOST_OR_DUMMY
                   自由表面粒子       → ..._IS_CONNECTED (探索の種)
                   内部粒子           → ..._IS_NOT_CONNECTED
              2. CONNECTED の粒子を見つけたら、その影響半径内にある
                 NOT_CONNECTED の粒子を CONNECTED に変える。
                 自分自身は CHECKED にして探索済みとする。
              3. 状態が変化しなくなるまで 2 を繰り返す。
              4. 最後まで NOT_CONNECTED のままの粒子があれば、
                 自由表面から孤立しているので警告を出す。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】FlagForCheckingBoundaryCondition[] を更新する
  【呼び出し元】exceptionalProcessingForBoundaryCondition()
  【計算量】最悪 O(N^2) × 繰り返し回数。かなり重い処理である。
=====================================================================*/
void checkBoundaryCondition( void ){
  int i,j,count;   /* i,j:粒子番号,
                      count:この周回で新たに探索した粒子の数。
                             0 になったら探索完了                    */
  double xij, yij, zij, distance2;  /* 相対位置と距離の 2 乗          */

  /* --- 手順 1 : 探索フラグの初期化 --- */
  for(i=0;i<NumberOfParticles;i++){
    if (BoundaryCondition[i]==GHOST_OR_DUMMY){
      FlagForCheckingBoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if (BoundaryCondition[i]==SURFACE_PARTICLE){
      /* 自由表面粒子そのものが探索の出発点になる                    */
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CONNECTED;
    }else{
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
    }
  }

  /* --- 手順 2, 3 : 変化がなくなるまで伝播を繰り返す --- */
  do {
    count=0;
    for(i=0;i<NumberOfParticles;i++){
      if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_CONNECTED){
	for(j=0;j<NumberOfParticles;j++){
	  if( j==i ) continue;
	  if((ParticleType[j]==GHOST) || (ParticleType[j]== DUMMY_WALL)) continue;

	  /* まだつながっていない粒子だけを調べる                    */
	  if(FlagForCheckingBoundaryCondition[j]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
	    xij = Position[j*3  ] - Position[i*3  ];
	    yij = Position[j*3+1] - Position[i*3+1];
	    zij = Position[j*3+2] - Position[i*3+2];
	    distance2 = (xij*xij)+(yij*yij)+(zij*zij);

	    /* 影響半径の外なら伝播しない
	       (2 乗同士の比較なので sqrt が不要で速い)              */
	    if(distance2>=Re2_forLaplacian)continue;

	    /* 影響半径内なので「つながった」ことにする              */
	    FlagForCheckingBoundaryCondition[j]=DIRICHLET_BOUNDARY_IS_CONNECTED;
	  }
	}
	/* 粒子 i からの伝播は済んだので探索済みにする               */
	FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CHECKED;
	count++;
      }
    }
  } while (count!=0); /* This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/
  /* count が 0 になる = 新たにつながった粒子がもう無い = 探索完了   */

  /* --- 手順 4 : 孤立した粒子があれば警告を出す --- */
  for(i=0;i<NumberOfParticles;i++){
    if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
      fprintf(stderr,"WARNING: There is no dirichlet boundary condition for %d-th particle.\n",i );
    }
  }
}


/*=====================================================================
  【関数名】increaseDiagonalTerm
  【機能】  自由表面 (ディリクレ境界) とつながっていない粒子について、
            係数行列の対角項を 2 倍にする。
            これにより行列が特異になるのを避け、解が発散しないように
            する救済処置である。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】CoefficientMatrix[] の対角項を書き換える
  【呼び出し元】exceptionalProcessingForBoundaryCondition()
=====================================================================*/
void increaseDiagonalTerm( void ){
  int i;                      /* 粒子番号                            */
  int n = NumberOfParticles;  /* 行列の一辺の大きさ                  */

  for(i=0;i<n;i++) {
    if(FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED ){
      CoefficientMatrix[i*n+i] = 2.0 * CoefficientMatrix[i*n+i];
    }
  }
}


/*=====================================================================
  【関数名】solveSimultaniousEquationsByGaussEliminationMethod
  【機能】  連立一次方程式 [A]{P} = {b} をガウスの消去法で解き、
            各粒子の圧力 Pressure[] を求める。

            (a) 前進消去 : 行 i を使って、その下の行 j (j > i) の
                i 列目を 0 にする。同じ操作を右辺 SourceTerm[] にも
                施して係数行列を上三角行列に変形する。
            (b) 後退代入 : 一番下の行から順に未知数を求める。
                  Pi = (bi − Σ_{j>i} A[i][j]・Pj) / A[i][i]

            内部粒子以外の行は飛ばしているため、
            自由表面粒子とダミー粒子の圧力は 0.0 のままになる。
            これが P = 0 のディリクレ境界条件に相当する。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Pressure[] を求め、CoefficientMatrix[] と SourceTerm[] を
            消去の過程で破壊的に書き換える
  【呼び出し元】calPressure()
  【計算量】O(N^3)。このプログラムで最も時間がかかる部分である。
            実用的な規模では ICCG 法などの反復解法に置き換える。
  【注意】  ピボット選択 (行の入れ替え) を行っていない。
            対角項に圧縮率項を足して対角優位にしてあるため、
            本プログラムの条件では問題なく解ける。
=====================================================================*/
void solveSimultaniousEquationsByGaussEliminationMethod( void ){
  int    i,j,k;       /* i:基準行, j:消去される行, k:列              */
  double c;           /* 消去の倍率 c = A[j][i] / A[i][i]            */
  double sumOfTerms;  /* 後退代入で使う Σ A[i][j]・Pj の値          */
  int    n = NumberOfParticles;  /* 未知数 (= 粒子) の個数           */

  /* 圧力を 0 で初期化する。
     内部粒子以外はこの 0 が最終的な圧力になる                       */
  for(i=0; i<n; i++){
    Pressure[i] = 0.0;
  }

  /* --- (a) 前進消去 --- */
  for(i=0; i<n-1; i++){
    /* 内部粒子の行だけを基準行として使う                            */
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;

    for(j=i+1; j<n; j++){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;

      /* 行 j の i 列目を 0 にするための倍率を求める                 */
      c = CoefficientMatrix[j*n+i]/CoefficientMatrix[i*n+i];

      /* 行 j から 行 i の c 倍を引く (i 列より右だけ計算すればよい) */
      for(k=i+1; k<n; k++){
	CoefficientMatrix[j*n+k] -= c * CoefficientMatrix[i*n+k];
      }
      /* 右辺にも同じ操作を施す                                      */
      SourceTerm[j] -= c*SourceTerm[i];
    }
  }

  /* --- (b) 後退代入 : 下の行から順に未知数を確定させる --- */
  for( i=n-1; i>=0; i--){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;

    sumOfTerms = 0.0;
    for( j=i+1; j<n; j++ ){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      /* すでに求まっている Pj を使って既知項をまとめる              */
      sumOfTerms += CoefficientMatrix[i*n+j] * Pressure[j];
    }
    Pressure[i] = (SourceTerm[i] - sumOfTerms)/CoefficientMatrix[i*n+i];
  }
}


/*=====================================================================
  【関数名】removeNegativePressure
  【機能】  負の圧力を 0 にクリップする。

            負圧をそのまま使うと、粒子どうしが引き合う力 (引張力) が
            働いて粒子が凝集し、自由表面が不自然に崩れてしまう。
            水は引張力に耐えられないという物理的事実にも合うため、
            負圧を切り捨てる。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Pressure[] を書き換える
  【呼び出し元】calPressure()
=====================================================================*/
void removeNegativePressure( void ){
  int i;   /* 粒子番号                                              */

  for(i=0;i<NumberOfParticles;i++) {
    if(Pressure[i]<0.0)Pressure[i]=0.0;
  }
}


/*=====================================================================
  【関数名】setMinimumPressure
  【機能】  各粒子 i について、自分自身と影響半径内の近傍粒子の中で
            最も小さい圧力を求め、MinimumPressure[i] に記録する。

            圧力勾配モデルで Pj − Pi の代わりに Pj − Pmin,i を使う
            ためのもの。こうすると (Pj − Pmin,i) が常に 0 以上と
            なるので、粒子間には必ず斥力 (反発力) だけが働き、
            引力による粒子の凝集や計算の破綻を防げる。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】MinimumPressure[] を更新する
  【呼び出し元】calPressure()  (圧力計算の最後)
  【計算量】O(N^2)
=====================================================================*/
void setMinimumPressure( void ){
  double xij, yij, zij, distance2;  /* 相対位置と距離の 2 乗          */
  int i,j;                          /* i:着目粒子, j:近傍粒子         */

  for(i=0;i<NumberOfParticles;i++) {
    /* 圧力を持たない粒子は対象外                                    */
    if(ParticleType[i]==GHOST || ParticleType[i]==DUMMY_WALL)continue;

    /* まず自分自身の圧力を最小値の初期値とする                      */
    MinimumPressure[i]=Pressure[i];

    for(j=0;j<NumberOfParticles;j++) {
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      /* ダミー壁は圧力を解いていない (常に 0) ので除外する。
	 含めると最小値が必ず 0 になってしまう                      */
      if(ParticleType[j]==DUMMY_WALL) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);

      /* 勾配モデル用の影響半径の外なら対象外                        */
      if(distance2>=Re2_forGradient)continue;

      /* より小さい圧力が見つかったら更新する                        */
      if( MinimumPressure[i] > Pressure[j] ){
	MinimumPressure[i] = Pressure[j];
      }
    }
  }
}


/*=====================================================================
  【関数名】calPressureGradient
  【機能】  圧力勾配 ∇P を MPS 法の勾配モデルで離散化し、
            そこから生じる加速度を求める。

            勾配モデル :
              <∇P>i = (DIM/n0) ・ Σ_j
                      (Pj − Pmin,i)/|rij|^2 ・ rij ・ w(rij)

            加速度は運動方程式  a = −(1/ρ)∇P  から求める。

            Pi の代わりに Pmin,i を使うのは setMinimumPressure() の
            説明のとおり、粒子間に斥力だけを働かせるためである。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Acceleration[] を「代入」で上書きする
            (moveParticle() でクリアされているので代入でよい)
  【呼び出し元】mainLoopOfSimulation()  (calPressure の直後)
  【計算量】O(N^2)
=====================================================================*/
void calPressureGradient( void ){
  int    i,j;    /* i:着目粒子, j:近傍粒子                           */
  double gradient_x, gradient_y, gradient_z;
                 /* 圧力勾配ベクトルの各成分 (総和を貯める変数)      */
  double xij, yij, zij;        /* 相対位置ベクトル rij               */
  double distance, distance2;  /* 粒子間距離とその 2 乗              */
  double w;                    /* 重み w(rij)                        */
  double pij;                  /* (Pj − Pmin,i)/|rij|^2              */
  double a;                    /* 勾配モデルの係数 a = DIM/n0        */

  a =DIM/N0_forGradient;

  for(i=0;i<NumberOfParticles;i++){
    /* 動かすのは流体粒子だけ                                        */
    if(ParticleType[i] != FLUID) continue;

    gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;

    for(j=0;j<NumberOfParticles;j++){
      if( j==i ) continue;
      if( ParticleType[j]==GHOST ) continue;
      /* ダミー壁は圧力を解いていないので勾配計算から除外する        */
      if( ParticleType[j]==DUMMY_WALL ) continue;

      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);

      if(distance<Re_forGradient){
	w =  weight(distance, Re_forGradient);
	/* 圧力差を距離の 2 乗で割る。
	   後で rij を掛けるので、全体として単位ベクトル × (差/距離)
	   という形になる                                            */
	pij = (Pressure[j] - MinimumPressure[i])/distance2;
	gradient_x += xij*pij*w;
	gradient_y += yij*pij*w;
	gradient_z += zij*pij*w;
      }
    }

    /* 総和に係数 DIM/n0 を掛けて圧力勾配にする                      */
    gradient_x *= a;
    gradient_y *= a;
    gradient_z *= a;

    /* 運動方程式 a = −(1/ρ)∇P より加速度を求める。
       ここは += ではなく = (代入) であることに注意。
       重力・粘性の加速度は moveParticle() で反映済みであり、
       ここでは圧力による加速度だけを改めて入れ直している           */
    Acceleration[i*3  ]= (-1.0)*gradient_x/FluidDensity;
    Acceleration[i*3+1]= (-1.0)*gradient_y/FluidDensity;
    Acceleration[i*3+2]= (-1.0)*gradient_z/FluidDensity;
  }
}


/*=====================================================================
  【関数名】moveParticleUsingPressureGradient
  【機能】  圧力勾配による加速度を使って、仮の速度・位置を修正し、
            次ステップの値として確定させる。

              u^(k+1) = u*  + a・Δt
              r^(k+1) = r*  + a・Δt^2

            位置の修正に u ではなく a・Δt^2 を使っているのは、
            「速度の修正量 (a・Δt) が Δt の間だけ作用した分だけ
              位置がずれる」と考えているためである。

            更新後は加速度を 0 にクリアし、次ステップに備える。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】Velocity[], Position[] を更新し、Acceleration[] を 0 にする
  【呼び出し元】mainLoopOfSimulation()  (1 ステップの最後)
=====================================================================*/
void moveParticleUsingPressureGradient( void ){
  int i;   /* 粒子番号                                              */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      /* 速度の修正 : u^(k+1) = u* + a・Δt                          */
      Velocity[i*3  ] +=Acceleration[i*3  ]*DT;
      Velocity[i*3+1] +=Acceleration[i*3+1]*DT;
      Velocity[i*3+2] +=Acceleration[i*3+2]*DT;

      /* 位置の修正 : r^(k+1) = r* + a・Δt^2                        */
      Position[i*3  ] +=Acceleration[i*3  ]*DT*DT;
      Position[i*3+1] +=Acceleration[i*3+1]*DT*DT;
      Position[i*3+2] +=Acceleration[i*3+2]*DT*DT;
    }
    /* 次ステップのために加速度をクリアする                          */
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
  }
}


/*=====================================================================
  【関数名】writeData_inProfFormat
  【機能】  計算結果をテキスト形式 (.prof) で出力する。
            自作の後処理プログラムや Excel などで読みやすい形式。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】ファイル "output_%04d.prof" を作成し、FileNumber を +1 する
  【呼び出し元】mainLoopOfSimulation()

  【出力ファイルの書式】
      1 行目        : 時刻 Time [s]
      2 行目        : 粒子数 NumberOfParticles
      3 行目以降    : 粒子 1 個につき 1 行、以下の 9 項目を空白区切り
                      粒子種類, 位置x, 位置y, 位置z,
                      速度x, 速度y, 速度z, 圧力, 粒子数密度
=====================================================================*/
void writeData_inProfFormat( void ){
  int i;                /* 粒子番号                                 */
  FILE *fp;             /* 出力ファイルのファイルポインタ            */
  char fileName[256];   /* 出力ファイル名を組み立てるバッファ        */

  /* 通し番号付きのファイル名を作る (例: output_0012.prof)          */
  sprintf(fileName, "output_%04d.prof",FileNumber);

  fp = fopen(fileName, "w");
  /* 【VC++対応の追加】ファイルが開けなかった場合の保護。
     チェックしないと、以降の fprintf で NULL ポインタを参照して
     異常終了してしまう                                             */
  if( fp == NULL ){
    fprintf(stderr,"ERROR: cannot open the file \"%s\".\n", fileName);
    exit(1);
  }

  fprintf(fp,"%lf\n",Time);                /* 1 行目: 時刻           */
  fprintf(fp,"%d\n",NumberOfParticles);    /* 2 行目: 粒子数         */

  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n"
	    ,ParticleType[i], Position[i*3], Position[i*3+1], Position[i*3+2]
	    ,Velocity[i*3], Velocity[i*3+1], Velocity[i*3+2], Pressure[i], NumberDensity[i]);
  }

  fclose(fp);
  FileNumber++;   /* 次回の出力に備えて番号を進める                  */
}


/*=====================================================================
  【関数名】writeData_inVtuFormat
  【機能】  計算結果を VTK の XML 形式 (.vtu = UnstructuredGrid) で
            出力する。ParaView などの可視化ソフトでそのまま開ける。

            各粒子を「1 点だけからなるセル (VTK_VERTEX, type=1)」と
            して登録し、点データとして以下を付ける。
              ParticleType : 粒子の種類 (整数)
              Velocity     : 速度の大きさ (スカラー)
              pressure     : 圧力

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】ファイル "particle_%04d.vtu" を作成する。
            ※ FileNumber の増加は writeData_inProfFormat() が行う。
              main ループでは vtu → prof の順に呼ぶことで、
              両者に同じ番号が振られるようになっている。
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void writeData_inVtuFormat( void ){
  int i;                          /* 粒子番号                       */
  double absoluteValueOfVelocity; /* 速度ベクトルの大きさ |u| [m/s] */
  FILE *fp;                       /* 出力ファイルのファイルポインタ  */
  char fileName[1024];            /* 出力ファイル名のバッファ        */

  /* 通し番号付きのファイル名を作る (例: particle_0012.vtu)          */
  sprintf(fileName, "particle_%04d.vtu", FileNumber);

  fp=fopen(fileName,"w");
  /* 【VC++対応の追加】ファイルが開けなかった場合の保護              */
  if( fp == NULL ){
    fprintf(stderr,"ERROR: cannot open the file \"%s\".\n", fileName);
    exit(1);
  }

  /* --- XML のヘッダ部 --- */
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  /* 粒子 1 個 = セル 1 個 = 点 1 個なので、どちらも粒子数と同じ      */
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",NumberOfParticles,NumberOfParticles);

  /* --- 座標 (Points) : x y z を粒子数分書き出す --- */
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",Position[i*3],Position[i*3+1],Position[i*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");

  /* --- 点データ (PointData) : 可視化で色を付ける量 --- */
  fprintf(fp,"<PointData>\n");

  /* 粒子の種類 (流体 / 壁 などを色分けするために使う)               */
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",ParticleType[i]);
  }
  fprintf(fp,"</DataArray>\n");

  /* 速度の大きさ |u| = sqrt(ux^2 + uy^2 + uz^2)                     */
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( Velocity[i*3]*Velocity[i*3] + Velocity[i*3+1]*Velocity[i*3+1] + Velocity[i*3+2]*Velocity[i*3+2] );
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }
  fprintf(fp,"</DataArray>\n");

  /* 圧力                                                            */
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)Pressure[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  /* --- セル情報 (Cells) --- */
  fprintf(fp,"<Cells>\n");

  /* connectivity : 各セルが使う点の番号。1 点なので i そのもの      */
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");

  /* offsets : connectivity 配列の区切り位置。1 点ずつなので i+1     */
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i+1);
  }
  fprintf(fp,"</DataArray>\n");

  /* types : セルの種類。1 は VTK_VERTEX (点) を表す                 */
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  /* --- XML の終端 --- */
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");

  fclose(fp);
}
