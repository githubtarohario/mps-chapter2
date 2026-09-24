/*=====================================================================
  mps.c
  (c) Kazuya SHIBATA, Kohei MUROTANI and Seiichi KOSHIZUKA (2014)

   Fluid Simulation Program Based on a Particle Method (the MPS method)
   Last update: May 21, 2014
=======================================================================*/
/*=====================================================================
 【このファイルについて】
   mps.c のオリジナル版 (書籍付属のまま) のバックアップ。
   コードには手を加えず、日本語のコメントだけを追加している。

 【プログラムの概要】
   MPS 法 (Moving Particle Semi-implicit method / 半陰的粒子法) で
   ダムブレイク (ダム崩壊) 問題を解く。
   箱の左側に置いた水柱 (幅 0.25m × 高さ 0.50m) が重力で崩れ、
   右へ流れ広がる様子を FINISH_TIME 秒まで計算する。

 【1 タイムステップの流れ】  mainLoopOfSimulation() を参照
   (1) 陽的計算 : 重力 → 粘性 → 仮の移動 → 衝突処理
   (2) 陰的計算 : 圧力のポアソン方程式を解いて圧力を求める
   (3) 修正     : 圧力勾配で速度・位置を修正する

 【出力ファイル】
   output_%04d.prof   : テキスト形式 (時刻・粒子数・各粒子の状態量)
   particle_%04d.vtu  : VTK 形式 (ParaView などで可視化できる)
=======================================================================*/
#include <stdio.h>    /* printf, fprintf, fopen などの入出力 */
#include <stdlib.h>   /* 標準ライブラリ (本プログラムでは実質未使用) */
#include <math.h>     /* sqrt (粒子間距離の計算) */
#include <string.h>   /* 文字列処理 (本プログラムでは実質未使用) */

/*---------------------------------------------------------------------
  【計算条件の設定 (2 次元計算用)】
---------------------------------------------------------------------*/
#define DIM                  2      /* 次元数 (2: 2 次元計算, 3: 3 次元計算) */
#define PARTICLE_DISTANCE    0.025  /* 初期粒子間距離 l0 [m] */
#define DT                   0.001  /* 時間刻み幅 Δt [s] */
#define OUTPUT_INTERVAL      20     /* 何ステップごとにファイル出力するか */

/* for three-dimensional simulation */
/* 3 次元計算をするときは、上の 4 行の代わりに下の 4 行を有効にする */
/*
#define DIM                  3
#define PARTICLE_DISTANCE    0.075
#define DT                   0.003
#define OUTPUT_INTERVAL      2
*/

/*---------------------------------------------------------------------
  【物理定数・計算パラメータ】
---------------------------------------------------------------------*/
#define ARRAY_SIZE           5000        /* 粒子数の上限 (配列の大きさ) */
#define FINISH_TIME          2.0         /* 計算を終える時刻 [s] */
#define KINEMATIC_VISCOSITY  (1.0E-6)    /* 動粘性係数 ν [m^2/s] (水) */
#define FLUID_DENSITY        1000.0      /* 流体の密度 ρ [kg/m^3] (水) */
#define G_X  0.0                         /* 重力加速度の x 成分 [m/s^2] */
#define G_Y  -9.8                        /* 重力加速度の y 成分 [m/s^2] (下向き) */
#define G_Z  0.0                         /* 重力加速度の z 成分 [m/s^2] */
#define RADIUS_FOR_NUMBER_DENSITY  (2.1*PARTICLE_DISTANCE) /* 粒子数密度の影響半径 re [m] */
#define RADIUS_FOR_GRADIENT        (2.1*PARTICLE_DISTANCE) /* 勾配モデルの影響半径 re [m] */
#define RADIUS_FOR_LAPLACIAN       (3.1*PARTICLE_DISTANCE) /* ラプラシアンモデルの影響半径 re [m] */
#define COLLISION_DISTANCE         (0.5*PARTICLE_DISTANCE) /* これより近づいたら衝突とみなす距離 [m] */
#define THRESHOLD_RATIO_OF_NUMBER_DENSITY  0.97   /* 自由表面判定の閾値 β (n < β・n0 なら表面粒子) */
#define COEFFICIENT_OF_RESTITUTION 0.2   /* 衝突時の反発係数 e */
#define COMPRESSIBILITY (0.45E-9)        /* 圧縮率 [1/Pa] (行列の対角項を少し大きくして安定化) */
#define EPS             (0.01 * PARTICLE_DISTANCE)  /* 初期配置の領域判定で使う微小量 (丸め誤差よけ) */
#define ON              1                /* フラグ: オン */
#define OFF             0                /* フラグ: オフ */
#define RELAXATION_COEFFICIENT_FOR_PRESSURE 0.2  /* 圧力計算の緩和係数 γ (ソース項に掛ける) */

/* 粒子の種類 (ParticleType[] に入る値) */
#define GHOST  -1        /* 計算対象外の粒子 (本プログラムでは生成されない) */
#define FLUID   0        /* 流体粒子 (運動方程式を解く) */
#define WALL    2        /* 壁粒子 (動かないが圧力は解く) */
#define DUMMY_WALL  3    /* ダミー壁粒子 (壁の外側に置き、粒子数密度を補うだけ。圧力は解かない) */

/* 圧力計算での境界条件 (BoundaryCondition[] に入る値) */
#define GHOST_OR_DUMMY  -1    /* ゴーストまたはダミー壁 (圧力計算に参加しない) */
#define SURFACE_PARTICLE 1    /* 自由表面の粒子 (圧力 0 のディリクレ境界) */
#define INNER_PARTICLE   0    /* 内部の粒子 (ポアソン方程式を解く) */

/* ディリクレ境界の連結性チェック (FlagForCheckingBoundaryCondition[] に入る値) */
#define DIRICHLET_BOUNDARY_IS_NOT_CONNECTED 0   /* 表面粒子とつながっていない (まだ見つかっていない) */
#define DIRICHLET_BOUNDARY_IS_CONNECTED     1   /* 表面粒子とつながっている (これから近傍を調べる) */
#define DIRICHLET_BOUNDARY_IS_CHECKED       2   /* つながっていて、近傍も調べ終えた */

/*---------------------------------------------------------------------
  【関数プロトタイプ宣言】
---------------------------------------------------------------------*/
void initializeParticlePositionAndVelocity_for2dim( void );
void initializeParticlePositionAndVelocity_for3dim( void );
void calConstantParameter( void );
void calNZeroAndLambda( void );
double weight( double distance, double re );
void mainLoopOfSimulation( void );
void calGravity( void );
void calViscosity( void );
void moveParticle( void );
void collision( void );
void calPressure( void );
void calNumberDensity( void );
void setBoundaryCondition( void );
void setSourceTerm( void );
void setMatrix( void );
void exceptionalProcessingForBoundaryCondition( void );
void checkBoundaryCondition( void );
void increaseDiagonalTerm( void );
void solveSimultaniousEquationsByGaussEliminationMethod( void );
void removeNegativePressure( void );
void setMinimumPressure( void );
void calPressureGradient( void );
void moveParticleUsingPressureGradient( void );
void writeData_inProfFormat( void );
void writeData_inVtuFormat( void );

/*---------------------------------------------------------------------
  【グローバル変数 (粒子ごとの配列)】
   ベクトル量は 1 粒子あたり 3 要素 (x, y, z) を並べて格納する。
   粒子 i の x 成分は [i*3], y 成分は [i*3+1], z 成分は [i*3+2]。
---------------------------------------------------------------------*/
static double Acceleration[3*ARRAY_SIZE];   /* 加速度 (x,y,z) [m/s^2] */
static int    ParticleType[ARRAY_SIZE];     /* 粒子の種類 (FLUID / WALL / DUMMY_WALL / GHOST) */
static double Position[3*ARRAY_SIZE];       /* 位置 (x,y,z) [m] */
static double Velocity[3*ARRAY_SIZE];       /* 速度 (x,y,z) [m/s] */
static double Pressure[ARRAY_SIZE];         /* 圧力 [Pa] */
static double NumberDensity[ARRAY_SIZE];    /* 粒子数密度 n (近傍粒子の重みの合計) */
static int    BoundaryCondition[ARRAY_SIZE];/* 圧力計算の境界条件 (INNER / SURFACE / GHOST_OR_DUMMY) */
static double SourceTerm[ARRAY_SIZE];       /* ポアソン方程式の右辺ベクトル {b} */
static int    FlagForCheckingBoundaryCondition[ARRAY_SIZE]; /* ディリクレ境界とつながっているかのチェック用フラグ */
static double CoefficientMatrix[ARRAY_SIZE * ARRAY_SIZE];   /* ポアソン方程式の係数行列 [A] (n×n, 行優先) */
static double MinimumPressure[ARRAY_SIZE];  /* 各粒子の近傍で最小の圧力 (圧力勾配の計算に使う) [Pa] */

/*---------------------------------------------------------------------
  【グローバル変数 (計算全体で使う値)】
---------------------------------------------------------------------*/
int    FileNumber;          /* 出力ファイルの通し番号 (output_%04d / particle_%04d) */
double Time;                /* 現在のシミュレーション時刻 [s] */
int    NumberOfParticles;   /* 粒子の総数 (流体 + 壁 + ダミー壁) */
double Re_forNumberDensity,Re2_forNumberDensity; /* 粒子数密度の影響半径 re [m] とその 2 乗 */
double Re_forGradient,     Re2_forGradient;      /* 勾配モデルの影響半径 re [m] とその 2 乗 */
double Re_forLaplacian,    Re2_forLaplacian;     /* ラプラシアンモデルの影響半径 re [m] とその 2 乗 */
double N0_forNumberDensity; /* 基準粒子数密度 n0 (粒子数密度用の影響半径で計算) */
double N0_forGradient;      /* 基準粒子数密度 n0 (勾配モデル用の影響半径で計算) */
double N0_forLaplacian;     /* 基準粒子数密度 n0 (ラプラシアンモデル用の影響半径で計算) */
double Lambda;              /* ラプラシアンモデルの係数 λ (距離の 2 乗の重み付き平均) */
double collisionDistance,collisionDistance2;     /* 衝突判定距離 [m] とその 2 乗 */
double FluidDensity;        /* 流体の密度 ρ [kg/m^3] */


/*=====================================================================
  【関数名】main
  【機能】  プログラムの入口。粒子の初期配置 → 定数の計算 →
            時間発展のメインループ、の順に実行する。
  【引数】  int    argc : コマンドライン引数の個数 (未使用)
            char** argv : コマンドライン引数の文字列配列 (未使用)
  【戻り値】int : 正常終了なら 0
  【呼び出し元】C ランタイム (プログラム起動時)
=====================================================================*/
int main( int argc, char** argv ) {

  printf("\n*** START PARTICLE-SIMULATION ***\n");
  if( DIM == 2 ){
    initializeParticlePositionAndVelocity_for2dim();
  }else{
    initializeParticlePositionAndVelocity_for3dim();
  }
  calConstantParameter();
  mainLoopOfSimulation();
  printf("*** END ***\n\n");
  return 0;
}


/*=====================================================================
  【関数名】initializeParticlePositionAndVelocity_for2dim
  【機能】  2 次元ダムブレイク問題の初期粒子配置を作る。
            格子状に候補点を並べ、領域ごとに粒子の種類を決める。
            (ダミー壁 → 壁 → 空洞 → 流体 の順に上書きしていく)
            全粒子の速度は 0 にする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】ParticleType[], Position[], Velocity[],
                              NumberOfParticles
  【呼び出し元】main()  (DIM == 2 のとき)
=====================================================================*/
void initializeParticlePositionAndVelocity_for2dim( void ){
  int iX, iY;                    /* 格子点の番号 (x 方向, y 方向) */
  int nX, nY;                    /* 格子点の数 (x 方向, y 方向) */
  double x, y, z;                /* 格子点の座標 [m] */
  int i = 0;                     /* 生成した粒子の番号 (最後に粒子数になる) */
  int flagOfParticleGeneration;  /* この格子点に粒子を置くかどうか (ON / OFF) */

  nX = (int)(1.0/PARTICLE_DISTANCE)+5;
  nY = (int)(0.6/PARTICLE_DISTANCE)+5;
  for(iX= -4;iX<nX;iX++){
    for(iY= -4;iY<nY;iY++){
      x = PARTICLE_DISTANCE * (double)(iX);
      y = PARTICLE_DISTANCE * (double)(iY);
      z = 0.0;
      flagOfParticleGeneration = OFF;

      /* dummy wall region */
      /* ダミー壁の領域 (壁の外側 4 層分) */
      if( ((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-4.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)) ){
	ParticleType[i]=DUMMY_WALL;
	flagOfParticleGeneration = ON;
      }

      /* wall region */
      /* 壁の領域 (内側 2 層分) */
      if( ((x>-2.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+2.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)) ){
	ParticleType[i]=WALL;
	flagOfParticleGeneration = ON;
      }

      /* wall region */
      /* 壁の領域 (上端の 2 層分) */
      if( ((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.6-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)) ){
	ParticleType[i]=WALL;
	flagOfParticleGeneration = ON;
      }

      /* empty region */
      /* 箱の中 (空洞) には粒子を置かない */
      if( ((x>0.0+EPS)&&(x<=1.00+EPS))&&( y>0.0+EPS )){
	flagOfParticleGeneration = OFF;
      }

      /* fluid region */
      /* 水柱の領域 (幅 0.25m × 高さ 0.50m) */
      if( ((x>0.0+EPS)&&(x<=0.25+EPS)) &&((y>0.0+EPS)&&(y<=0.50+EPS)) ){
	ParticleType[i]=FLUID;
	flagOfParticleGeneration = ON;
      }

      if( flagOfParticleGeneration == ON){
	Position[i*3]=x; Position[i*3+1]=y; Position[i*3+2]=z;
	i++;
      }
    }
  }
  NumberOfParticles = i;
  for(i=0;i<NumberOfParticles*3;i++) { Velocity[i]=0.0; }
}


/*=====================================================================
  【関数名】initializeParticlePositionAndVelocity_for3dim
  【機能】  3 次元ダムブレイク問題の初期粒子配置を作る。
            2 次元版に z 方向 (奥行き 0.3m) を加えたもの。
            全粒子の速度は 0 にする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】ParticleType[], Position[], Velocity[],
                              NumberOfParticles
  【呼び出し元】main()  (DIM != 2 のとき)
=====================================================================*/
void initializeParticlePositionAndVelocity_for3dim( void ){
  int iX, iY, iZ;                /* 格子点の番号 (x, y, z 方向) */
  int nX, nY, nZ;                /* 格子点の数 (x, y, z 方向) */
  double x, y, z;                /* 格子点の座標 [m] */
  int i = 0;                     /* 生成した粒子の番号 (最後に粒子数になる) */
  int flagOfParticleGeneration;  /* この格子点に粒子を置くかどうか (ON / OFF) */

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
	/* ダミー壁の領域 (壁の外側 4 層分) */
	if( (((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-4.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-4.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+4.0*PARTICLE_DISTANCE+EPS ))){
	  ParticleType[i]=DUMMY_WALL;
	  flagOfParticleGeneration = ON;
	}

	/* wall region */
	/* 壁の領域 (内側 2 層分) */
	if( (((x>-2.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+2.0*PARTICLE_DISTANCE+EPS))&&( (y>0.0-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-2.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+2.0*PARTICLE_DISTANCE+EPS ))){
	  ParticleType[i]=WALL;
	  flagOfParticleGeneration = ON;
	}

	/* wall region */
	/* 壁の領域 (上端の 2 層分) */
	if( (((x>-4.0*PARTICLE_DISTANCE+EPS)&&(x<=1.00+4.0*PARTICLE_DISTANCE+EPS))&&( (y>0.6-2.0*PARTICLE_DISTANCE+EPS )&&(y<=0.6+EPS)))&&( (z>0.0-4.0*PARTICLE_DISTANCE+EPS)&&(z<=0.3+4.0*PARTICLE_DISTANCE+EPS ))){
	  ParticleType[i]=WALL;
	  flagOfParticleGeneration = ON;
	}

	/* empty region */
	/* 箱の中 (空洞) には粒子を置かない */
	if( (((x>0.0+EPS)&&(x<=1.00+EPS))&&( y>0.0+EPS ))&&( (z>0.0+EPS )&&(z<=0.3+EPS ))){
	  flagOfParticleGeneration = OFF;
	}

	/* fluid region */
	/* 水柱の領域 (幅 0.25m × 高さ 0.50m × 奥行き 0.3m) */
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
  NumberOfParticles = i;
  for(i=0;i<NumberOfParticles*3;i++) { Velocity[i]=0.0; }
}


/*=====================================================================
  【関数名】calConstantParameter
  【機能】  計算中に変わらない定数 (影響半径・基準粒子数密度・λ・
            密度・衝突距離など) を計算し、時刻と出力番号を 0 にする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Re_for*, Re2_for*, N0_for*, Lambda,
            FluidDensity, collisionDistance, collisionDistance2,
            FileNumber, Time
  【呼び出し元】main()
=====================================================================*/
void calConstantParameter( void ){

  Re_forNumberDensity  = RADIUS_FOR_NUMBER_DENSITY;
  Re_forGradient       = RADIUS_FOR_GRADIENT;
  Re_forLaplacian      = RADIUS_FOR_LAPLACIAN;
  Re2_forNumberDensity = Re_forNumberDensity*Re_forNumberDensity;
  Re2_forGradient      = Re_forGradient*Re_forGradient;
  Re2_forLaplacian     = Re_forLaplacian*Re_forLaplacian;
  calNZeroAndLambda();
  FluidDensity       = FLUID_DENSITY;
  collisionDistance  = COLLISION_DISTANCE;
  collisionDistance2 = collisionDistance*collisionDistance;
  FileNumber=0;
  Time=0.0;
}


/*=====================================================================
  【関数名】calNZeroAndLambda
  【機能】  基準粒子数密度 n0 (3 種類の影響半径それぞれ) と、
            ラプラシアンモデルの係数 λ を求める。
            原点に置いた粒子のまわりに、粒子を格子状に理想配置した
            ときの重みの合計として計算する。
            λ = Σ(r^2 × w(r)) / Σ w(r)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】N0_forNumberDensity, N0_forGradient,
                              N0_forLaplacian, Lambda
  【呼び出し元】calConstantParameter()
=====================================================================*/
void calNZeroAndLambda( void ){
  int iX, iY, iZ;               /* 周囲の格子点の番号 (x, y, z 方向) */
  int iZ_start, iZ_end;         /* z 方向のループ範囲 (2 次元では 0 だけ) */
  double xj, yj, zj, distance, distance2;  /* 周囲の粒子 j の座標 [m]、粒子 i との距離 [m] とその 2 乗 */
  double xi, yi, zi;            /* 中心の粒子 i の座標 [m] (原点) */

  if( DIM == 2 ){
    iZ_start = 0; iZ_end = 1;
  }else{
    iZ_start = -4; iZ_end = 5;
  }

  N0_forNumberDensity = 0.0;
  N0_forGradient      = 0.0;
  N0_forLaplacian     = 0.0;
  Lambda              = 0.0;
  xi = 0.0;  yi = 0.0;  zi = 0.0;

  for(iX= -4;iX<5;iX++){
    for(iY= -4;iY<5;iY++){
      for(iZ= iZ_start;iZ<iZ_end;iZ++){
	if( ((iX==0)&&(iY==0)) && (iZ==0) )continue;   /* 自分自身は数えない */
	xj = PARTICLE_DISTANCE * (double)(iX);
	yj = PARTICLE_DISTANCE * (double)(iY);
	zj = PARTICLE_DISTANCE * (double)(iZ);
	distance2 = (xj-xi)*(xj-xi)+(yj-yi)*(yj-yi)+(zj-zi)*(zj-zi);
	distance = sqrt(distance2);
	N0_forNumberDensity += weight(distance, Re_forNumberDensity);
	N0_forGradient      += weight(distance, Re_forGradient);
	N0_forLaplacian     += weight(distance, Re_forLaplacian);
	Lambda              += distance2 * weight(distance, Re_forLaplacian);
      }
    }
  }
  Lambda = Lambda/N0_forLaplacian;
}


/*=====================================================================
  【関数名】weight
  【機能】  MPS 法の重み関数 w(r) を計算する。
                w(r) = re/r - 1   (r <  re)
                w(r) = 0          (r >= re)
            近い粒子ほど重みが大きく、影響半径 re より遠い粒子は 0。
  【引数】  double distance : 2 粒子間の距離 r [m] (0 より大きいこと)
            double re       : 影響半径 re [m]
  【戻り値】double : 重み w(r)
  【呼び出し元】calNZeroAndLambda(), calViscosity(), calNumberDensity(),
                setMatrix(), calPressureGradient()
=====================================================================*/
double weight( double distance, double re ){
  double weightIJ;   /* 粒子 i と j の間の重み (戻り値) */

  if( distance >= re ){
    weightIJ = 0.0;
  }else{
    weightIJ = (re/distance) - 1.0;
  }
  return weightIJ;
}


/*=====================================================================
  【関数名】mainLoopOfSimulation
  【機能】  時間発展のメインループ。FINISH_TIME に達するまで
            1 ステップずつ計算を進め、OUTPUT_INTERVAL ステップごとに
            経過を表示してファイルに書き出す。
            (最初に初期状態も書き出す)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Time
  【呼び出し元】main()
=====================================================================*/
void mainLoopOfSimulation( void ){
  int iTimeStep = 0;   /* 現在のタイムステップ番号 */

  writeData_inVtuFormat();
  writeData_inProfFormat();

  while(1){
    calGravity();                          /* (1) 陽的計算: 重力        */
    calViscosity();                        /*               粘性        */
    moveParticle();                        /*               仮の移動    */
    collision();                           /*               衝突処理    */
    calPressure();                         /* (2) 陰的計算: 圧力        */
    calPressureGradient();                 /* (3) 修正:     圧力勾配    */
    moveParticleUsingPressureGradient();   /*               速度・位置の修正 */
    iTimeStep++;
    Time += DT;
    if( (iTimeStep % OUTPUT_INTERVAL) == 0 ){
      printf("TimeStepNumber: %4d   Time: %lf(s)   NumberOfParticless: %d\n", iTimeStep, Time, NumberOfParticles);
      writeData_inVtuFormat();
      writeData_inProfFormat();
    }
    if( Time >= FINISH_TIME ){break;}
  }
}


/*=====================================================================
  【関数名】calGravity
  【機能】  流体粒子の加速度に重力加速度 (G_X, G_Y, G_Z) を設定する。
            流体以外の粒子の加速度は 0 にする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Acceleration[]
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void calGravity( void ){
  int i;   /* 粒子番号 */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Acceleration[i*3  ]=G_X;
      Acceleration[i*3+1]=G_Y;
      Acceleration[i*3+2]=G_Z;
    }else{
      Acceleration[i*3  ]=0.0;
      Acceleration[i*3+1]=0.0;
      Acceleration[i*3+2]=0.0;
    }
  }
}


/*=====================================================================
  【関数名】calViscosity
  【機能】  粘性項 ν∇^2u を MPS 法のラプラシアンモデルで計算し、
            流体粒子の加速度に加える。
                ν∇^2u_i = ν × 2d/(λ n0) × Σ (u_j - u_i) w(r_ij)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Acceleration[]
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void calViscosity( void ){
  int i,j;                /* i: 注目する粒子の番号, j: 近傍粒子の番号 */
  double viscosityTerm_x, viscosityTerm_y, viscosityTerm_z;  /* 粒子 i の粘性項 (x,y,z) [m/s^2] */
  double distance, distance2;  /* 粒子 i と j の距離 [m] とその 2 乗 */
  double w;               /* 粒子 i と j の間の重み */
  double xij, yij, zij;   /* 粒子 i から j への相対位置 (x,y,z) [m] */
  double a;               /* ラプラシアンモデルの係数 ν×2d/(λ n0) */

  a = (KINEMATIC_VISCOSITY)*(2.0*DIM)/(N0_forLaplacian*Lambda);
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] != FLUID) continue;
    viscosityTerm_x = 0.0;  viscosityTerm_y = 0.0;  viscosityTerm_z = 0.0;

    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      if(distance<Re_forLaplacian){
	w =  weight(distance, Re_forLaplacian);
	viscosityTerm_x +=(Velocity[j*3  ]-Velocity[i*3  ])*w;
	viscosityTerm_y +=(Velocity[j*3+1]-Velocity[i*3+1])*w;
	viscosityTerm_z +=(Velocity[j*3+2]-Velocity[i*3+2])*w;
      }
    }
    viscosityTerm_x = viscosityTerm_x * a;
    viscosityTerm_y = viscosityTerm_y * a;
    viscosityTerm_z = viscosityTerm_z * a;
    Acceleration[i*3  ] += viscosityTerm_x;
    Acceleration[i*3+1] += viscosityTerm_y;
    Acceleration[i*3+2] += viscosityTerm_z;
  }
}


/*=====================================================================
  【関数名】moveParticle
  【機能】  重力と粘性による加速度で、流体粒子の速度と位置を仮に
            進める (仮の速度 u*, 仮の位置 r*)。
            使い終わった加速度は全粒子について 0 に戻す。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Velocity[], Position[], Acceleration[]
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void moveParticle( void ){
  int i;   /* 粒子番号 */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Velocity[i*3  ] += Acceleration[i*3  ]*DT;
      Velocity[i*3+1] += Acceleration[i*3+1]*DT;
      Velocity[i*3+2] += Acceleration[i*3+2]*DT;

      Position[i*3  ] += Velocity[i*3  ]*DT;
      Position[i*3+1] += Velocity[i*3+1]*DT;
      Position[i*3+2] += Velocity[i*3+2]*DT;
    }
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
  }
}


/*=====================================================================
  【関数名】collision
  【機能】  流体粒子がほかの粒子に COLLISION_DISTANCE より近づき、
            さらに互いに近づく向きに動いているとき、反発係数 e の
            衝突として速度を修正する (粒子どうしのめり込み防止)。
            速度の変化分だけ位置も修正する。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Velocity[], Position[]
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void collision( void ){
  int    i,j;                /* i: 注目する粒子の番号, j: 相手の粒子の番号 */
  double xij, yij, zij;      /* 粒子 i から j への相対位置 (x,y,z) [m] */
  double distance,distance2; /* 粒子 i と j の距離 [m] とその 2 乗 */
  double forceDT; /* forceDT is the impulse of collision between particles */
                  /* 衝突の力積 (最初は近づく向きの相対速度として使う) */
  double mi, mj;             /* 粒子 i, j の質量 (ここでは密度で代用) */
  double velocity_ix, velocity_iy, velocity_iz;  /* 衝突後の粒子 i の速度 (x,y,z) [m/s] */
  double e = COEFFICIENT_OF_RESTITUTION;         /* 反発係数 */
  static double VelocityAfterCollision[3*ARRAY_SIZE];  /* 全粒子の衝突後の速度 (x,y,z) [m/s] */

  for(i=0;i<3*NumberOfParticles;i++){
    VelocityAfterCollision[i] = Velocity[i];
  }
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      mi = FluidDensity;
      velocity_ix = Velocity[i*3  ];
      velocity_iy = Velocity[i*3+1];
      velocity_iz = Velocity[i*3+2];
      for(j=0;j<NumberOfParticles;j++){
	if( (j==i) || (ParticleType[j]==GHOST) ) continue;
	xij = Position[j*3  ] - Position[i*3  ];
	yij = Position[j*3+1] - Position[i*3+1];
	zij = Position[j*3+2] - Position[i*3+2];
	distance2 = (xij*xij) + (yij*yij) + (zij*zij);
	if(distance2<collisionDistance2){
	  distance = sqrt(distance2);
	  forceDT = (velocity_ix-Velocity[j*3  ])*(xij/distance)
	           +(velocity_iy-Velocity[j*3+1])*(yij/distance)
	           +(velocity_iz-Velocity[j*3+2])*(zij/distance);
	  if(forceDT > 0.0){   /* 互いに近づいているときだけ反発させる */
	    mj = FluidDensity;
	    forceDT *= (1.0+e)*mi*mj/(mi+mj);
	    velocity_ix -= (forceDT/mi)*(xij/distance);
	    velocity_iy -= (forceDT/mi)*(yij/distance);
	    velocity_iz -= (forceDT/mi)*(zij/distance);
	    /*
	    if(j>i){ fprintf(stderr,"WARNING: Collision occured between %d and %d particles.\n",i,j); }
	    */
	  }
	}
      }
      VelocityAfterCollision[i*3  ] = velocity_ix;
      VelocityAfterCollision[i*3+1] = velocity_iy;
      VelocityAfterCollision[i*3+2] = velocity_iz;
    }
  }
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Position[i*3  ] += (VelocityAfterCollision[i*3  ]-Velocity[i*3  ])*DT;
      Position[i*3+1] += (VelocityAfterCollision[i*3+1]-Velocity[i*3+1])*DT;
      Position[i*3+2] += (VelocityAfterCollision[i*3+2]-Velocity[i*3+2])*DT;
      Velocity[i*3  ] = VelocityAfterCollision[i*3  ];
      Velocity[i*3+1] = VelocityAfterCollision[i*3+1];
      Velocity[i*3+2] = VelocityAfterCollision[i*3+2];
    }
  }
}


/*=====================================================================
  【関数名】calPressure
  【機能】  圧力のポアソン方程式 [A]{P} = {b} を組み立てて解き、
            各粒子の圧力を求める。処理の順番を並べただけの関数。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void calPressure( void ){
  calNumberDensity();                                    /* 粒子数密度 n を数える */
  setBoundaryCondition();                                /* 自由表面か内部かを判定 */
  setSourceTerm();                                       /* 右辺ベクトル {b} を作る */
  setMatrix();                                           /* 係数行列 [A] を作る */
  solveSimultaniousEquationsByGaussEliminationMethod();  /* 連立方程式を解く */
  removeNegativePressure();                              /* 負の圧力を 0 にする */
  setMinimumPressure();                                  /* 近傍の最小圧力を記録 */
}


/*=====================================================================
  【関数名】calNumberDensity
  【機能】  各粒子の粒子数密度 n_i = Σ w(r_ij) を計算する。
            (影響半径は RADIUS_FOR_NUMBER_DENSITY)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】NumberDensity[]
  【呼び出し元】calPressure()
=====================================================================*/
void calNumberDensity( void ){
  int    i,j;                /* i: 注目する粒子の番号, j: 近傍粒子の番号 */
  double xij, yij, zij;      /* 粒子 i から j への相対位置 (x,y,z) [m] */
  double distance, distance2;/* 粒子 i と j の距離 [m] とその 2 乗 */
  double w;                  /* 粒子 i と j の間の重み */

  for(i=0;i<NumberOfParticles;i++){
    NumberDensity[i] = 0.0;
    if(ParticleType[i] == GHOST) continue;
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      w =  weight(distance, Re_forNumberDensity);
      NumberDensity[i] += w;
    }
  }
}


/*=====================================================================
  【関数名】setBoundaryCondition
  【機能】  各粒子の圧力計算での扱いを決める。
              ゴースト・ダミー壁       → GHOST_OR_DUMMY (計算しない)
              n_i < β × n0 の粒子     → SURFACE_PARTICLE (自由表面, 圧力 0)
              それ以外                 → INNER_PARTICLE (方程式を解く)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】BoundaryCondition[]
  【呼び出し元】calPressure()
=====================================================================*/
void setBoundaryCondition( void ){
  int i;                                          /* 粒子番号 */
  double n0 = N0_forNumberDensity;                /* 基準粒子数密度 n0 */
  double beta = THRESHOLD_RATIO_OF_NUMBER_DENSITY;/* 自由表面判定の閾値 β */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ){
      BoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if( NumberDensity[i] < beta * n0 ){
      BoundaryCondition[i]=SURFACE_PARTICLE;
    }else{
      BoundaryCondition[i]=INNER_PARTICLE;
    }
  }
}


/*=====================================================================
  【関数名】setSourceTerm
  【機能】  ポアソン方程式の右辺ベクトル {b} を作る。
            内部粒子:  b_i = γ × (1/Δt^2) × (n_i - n0)/n0
            表面粒子:  b_i = 0 (圧力 0 のディリクレ境界)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】SourceTerm[]
  【呼び出し元】calPressure()
=====================================================================*/
void setSourceTerm( void ){
  int i;                                              /* 粒子番号 */
  double n0    = N0_forNumberDensity;                 /* 基準粒子数密度 n0 */
  double gamma = RELAXATION_COEFFICIENT_FOR_PRESSURE; /* 緩和係数 γ */

  for(i=0;i<NumberOfParticles;i++){
    SourceTerm[i]=0.0;
    if(ParticleType[i]==GHOST || ParticleType[i]== DUMMY_WALL ) continue;
    if(BoundaryCondition[i]==INNER_PARTICLE){
      SourceTerm[i] = gamma * (1.0/(DT*DT))*((NumberDensity[i]-n0)/n0);
    }else if(BoundaryCondition[i]==SURFACE_PARTICLE){
      SourceTerm[i]=0.0;
    }
  }
}


/*=====================================================================
  【関数名】setMatrix
  【機能】  ポアソン方程式の係数行列 [A] を作る。
            内部粒子 i の行に、ラプラシアンモデルの係数
            2d/(λ n0) × w(r_ij) / ρ を設定する。
            非対角 A_ij = -係数、対角 A_ii = Σ係数 + 圧縮率/Δt^2。
            最後に、ディリクレ境界とつながらない粒子の例外処理を行う。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】CoefficientMatrix[]
  【呼び出し元】calPressure()
=====================================================================*/
void setMatrix( void ){
  double xij, yij, zij;        /* 粒子 i から j への相対位置 (x,y,z) [m] */
  double distance, distance2;  /* 粒子 i と j の距離 [m] とその 2 乗 */
  double coefficientIJ;        /* 粒子 i と j の間の係数 */
  double n0 = N0_forLaplacian; /* 基準粒子数密度 n0 (ラプラシアン用) */
  int    i,j;                  /* i: 行 (注目する粒子), j: 列 (近傍粒子) */
  double a;                    /* ラプラシアンモデルの係数 2d/(λ n0) */
  int n = NumberOfParticles;   /* 行列の大きさ (粒子数) */

  for(i=0;i<NumberOfParticles;i++){
    for(j=0;j<NumberOfParticles;j++){
      CoefficientMatrix[i*n+j] = 0.0;
    }
  }

  a = 2.0*DIM/(n0*Lambda);
  for(i=0;i<NumberOfParticles;i++){
    if(BoundaryCondition[i] != INNER_PARTICLE) continue;
    for(j=0;j<NumberOfParticles;j++){
      if( (j==i) || (BoundaryCondition[j]==GHOST_OR_DUMMY) ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);
      distance  = sqrt(distance2);
      if(distance>=Re_forLaplacian)continue;
      coefficientIJ = a * weight(distance, Re_forLaplacian)/FluidDensity;
      CoefficientMatrix[i*n+j]  = (-1.0)*coefficientIJ;
      CoefficientMatrix[i*n+i] += coefficientIJ;
    }
    CoefficientMatrix[i*n+i] += (COMPRESSIBILITY)/(DT*DT);
  }
  exceptionalProcessingForBoundaryCondition();
}


/*=====================================================================
  【関数名】exceptionalProcessingForBoundaryCondition
  【機能】  ディリクレ境界 (自由表面) とつながっていない粒子の集まり
            があると行列が解けなくなるため、その粒子を探して
            行列の対角項を大きくする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】setMatrix()
=====================================================================*/
void exceptionalProcessingForBoundaryCondition( void ){
  /* If tere is no Dirichlet boundary condition on the fluid,
     increase the diagonal terms of the matrix for an exception. This allows us to solve the matrix without Dirichlet boundary conditions. */
  checkBoundaryCondition();
  increaseDiagonalTerm();
}


/*=====================================================================
  【関数名】checkBoundaryCondition
  【機能】  自由表面粒子から出発し、影響半径 (ラプラシアン用) 内で
            つながっている粒子を次々にたどって印を付ける。
            最後まで印が付かなかった粒子は、ディリクレ境界と
            つながっていないので警告を表示する。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】FlagForCheckingBoundaryCondition[]
  【呼び出し元】exceptionalProcessingForBoundaryCondition()
=====================================================================*/
void checkBoundaryCondition( void ){
  int i,j,count;                     /* i, j: 粒子番号, count: この周回で調べ終えた粒子の数 */
  double xij, yij, zij, distance2;   /* 粒子 i から j への相対位置 (x,y,z) [m] と距離の 2 乗 */

  for(i=0;i<NumberOfParticles;i++){
    if (BoundaryCondition[i]==GHOST_OR_DUMMY){
      FlagForCheckingBoundaryCondition[i]=GHOST_OR_DUMMY;
    }else if (BoundaryCondition[i]==SURFACE_PARTICLE){
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CONNECTED;
    }else{
      FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_NOT_CONNECTED;
    }
  }

  do {
    count=0;
    for(i=0;i<NumberOfParticles;i++){
      if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_CONNECTED){
	for(j=0;j<NumberOfParticles;j++){
	  if( j==i ) continue;
	  if((ParticleType[j]==GHOST) || (ParticleType[j]== DUMMY_WALL)) continue;
	  if(FlagForCheckingBoundaryCondition[j]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
	    xij = Position[j*3  ] - Position[i*3  ];
	    yij = Position[j*3+1] - Position[i*3+1];
	    zij = Position[j*3+2] - Position[i*3+2];
	    distance2 = (xij*xij)+(yij*yij)+(zij*zij);
	    if(distance2>=Re2_forLaplacian)continue;
	    FlagForCheckingBoundaryCondition[j]=DIRICHLET_BOUNDARY_IS_CONNECTED;
	  }
	}
	FlagForCheckingBoundaryCondition[i]=DIRICHLET_BOUNDARY_IS_CHECKED;
	count++;
      }
    }
  } while (count!=0); /* This procedure is repeated until the all fluid or wall particles (which have Dirhchlet boundary condition in the particle group) are in the state of "DIRICHLET_BOUNDARY_IS_CHECKED".*/
                      /* 新しくつながった粒子がなくなるまで繰り返す */

  for(i=0;i<NumberOfParticles;i++){
    if(FlagForCheckingBoundaryCondition[i]==DIRICHLET_BOUNDARY_IS_NOT_CONNECTED){
      fprintf(stderr,"WARNING: There is no dirichlet boundary condition for %d-th particle.\n",i );
    }
  }
}


/*=====================================================================
  【関数名】increaseDiagonalTerm
  【機能】  ディリクレ境界とつながっていない粒子について、係数行列の
            対角項を 2 倍にする。行列が特異 (解けない) になるのを防ぐ。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】CoefficientMatrix[]
  【呼び出し元】exceptionalProcessingForBoundaryCondition()
=====================================================================*/
void increaseDiagonalTerm( void ){
  int i;                       /* 粒子番号 (行列の行) */
  int n = NumberOfParticles;   /* 行列の大きさ (粒子数) */

  for(i=0;i<n;i++) {
    if(FlagForCheckingBoundaryCondition[i] == DIRICHLET_BOUNDARY_IS_NOT_CONNECTED ){
      CoefficientMatrix[i*n+i] = 2.0 * CoefficientMatrix[i*n+i];
    }
  }
}


/*=====================================================================
  【関数名】solveSimultaniousEquationsByGaussEliminationMethod
  【機能】  連立一次方程式 [A]{P} = {b} をガウスの消去法
            (前進消去 → 後退代入) で解き、圧力を求める。
            内部粒子以外の圧力は 0 のまま。
            ※前進消去で CoefficientMatrix[] と SourceTerm[] は
              書き換えられる。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Pressure[], CoefficientMatrix[], SourceTerm[]
  【呼び出し元】calPressure()
=====================================================================*/
void solveSimultaniousEquationsByGaussEliminationMethod( void ){
  int    i,j,k;               /* i: ピボット行, j: 消去する行, k: 列 */
  double c;                   /* 前進消去で行 i に掛ける倍率 A_ji / A_ii */
  double sumOfTerms;          /* 後退代入での Σ A_ij × P_j */
  int    n = NumberOfParticles;  /* 行列の大きさ (粒子数) */

  for(i=0; i<n; i++){
    Pressure[i] = 0.0;
  }
  /* 前進消去 */
  for(i=0; i<n-1; i++){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;
    for(j=i+1; j<n; j++){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      c = CoefficientMatrix[j*n+i]/CoefficientMatrix[i*n+i];
      for(k=i+1; k<n; k++){
	CoefficientMatrix[j*n+k] -= c * CoefficientMatrix[i*n+k];
      }
      SourceTerm[j] -= c*SourceTerm[i];
    }
  }
  /* 後退代入 */
  for( i=n-1; i>=0; i--){
    if ( BoundaryCondition[i] != INNER_PARTICLE ) continue;
    sumOfTerms = 0.0;
    for( j=i+1; j<n; j++ ){
      if(BoundaryCondition[j]==GHOST_OR_DUMMY) continue;
      sumOfTerms += CoefficientMatrix[i*n+j] * Pressure[j];
    }
    Pressure[i] = (SourceTerm[i] - sumOfTerms)/CoefficientMatrix[i*n+i];
  }
}


/*=====================================================================
  【関数名】removeNegativePressure
  【機能】  負になった圧力を 0 にする (負圧による不安定を防ぐ)。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Pressure[]
  【呼び出し元】calPressure()
=====================================================================*/
void removeNegativePressure( void ){
  int i;   /* 粒子番号 */

  for(i=0;i<NumberOfParticles;i++) {
    if(Pressure[i]<0.0)Pressure[i]=0.0;
  }
}


/*=====================================================================
  【関数名】setMinimumPressure
  【機能】  各粒子について、自分と近傍粒子 (勾配用の影響半径内) の
            圧力の最小値を求める。圧力勾配の計算で、この最小値を
            基準にすることで粒子間に常に斥力が働くようにする。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】MinimumPressure[]
  【呼び出し元】calPressure()
=====================================================================*/
void setMinimumPressure( void ){
  double xij, yij, zij, distance2;  /* 粒子 i から j への相対位置 (x,y,z) [m] と距離の 2 乗 */
  int i,j;                          /* i: 注目する粒子の番号, j: 近傍粒子の番号 */

  for(i=0;i<NumberOfParticles;i++) {
    if(ParticleType[i]==GHOST || ParticleType[i]==DUMMY_WALL)continue;
    MinimumPressure[i]=Pressure[i];
    for(j=0;j<NumberOfParticles;j++) {
      if( (j==i) || (ParticleType[j]==GHOST) ) continue;
      if(ParticleType[j]==DUMMY_WALL) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij)+(yij*yij)+(zij*zij);
      if(distance2>=Re2_forGradient)continue;
      if( MinimumPressure[i] > Pressure[j] ){
	MinimumPressure[i] = Pressure[j];
      }
    }
  }
}


/*=====================================================================
  【関数名】calPressureGradient
  【機能】  流体粒子の圧力勾配を MPS 法の勾配モデルで計算し、
            加速度 -(1/ρ)∇P として設定する。
                ∇P_i = d/n0 × Σ (P_j - Pmin_i)/r_ij^2 × (r_j - r_i) w(r_ij)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Acceleration[]
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void calPressureGradient( void ){
  int    i,j;                  /* i: 注目する粒子の番号, j: 近傍粒子の番号 */
  double gradient_x, gradient_y, gradient_z;  /* 粒子 i の圧力勾配 (x,y,z) [Pa/m] */
  double xij, yij, zij;        /* 粒子 i から j への相対位置 (x,y,z) [m] */
  double distance, distance2;  /* 粒子 i と j の距離 [m] とその 2 乗 */
  double w,pij;                /* w: 重み, pij: (P_j - Pmin_i) / r_ij^2 */
  double a;                    /* 勾配モデルの係数 d/n0 */

  a =DIM/N0_forGradient;
  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] != FLUID) continue;
    gradient_x = 0.0;  gradient_y = 0.0;  gradient_z = 0.0;
    for(j=0;j<NumberOfParticles;j++){
      if( j==i ) continue;
      if( ParticleType[j]==GHOST ) continue;
      if( ParticleType[j]==DUMMY_WALL ) continue;
      xij = Position[j*3  ] - Position[i*3  ];
      yij = Position[j*3+1] - Position[i*3+1];
      zij = Position[j*3+2] - Position[i*3+2];
      distance2 = (xij*xij) + (yij*yij) + (zij*zij);
      distance = sqrt(distance2);
      if(distance<Re_forGradient){
	w =  weight(distance, Re_forGradient);
	pij = (Pressure[j] - MinimumPressure[i])/distance2;
	gradient_x += xij*pij*w;
	gradient_y += yij*pij*w;
	gradient_z += zij*pij*w;
      }
    }
    gradient_x *= a;
    gradient_y *= a;
    gradient_z *= a;
    Acceleration[i*3  ]= (-1.0)*gradient_x/FluidDensity;
    Acceleration[i*3+1]= (-1.0)*gradient_y/FluidDensity;
    Acceleration[i*3+2]= (-1.0)*gradient_z/FluidDensity;
  }
}


/*=====================================================================
  【関数名】moveParticleUsingPressureGradient
  【機能】  圧力勾配による加速度で、流体粒子の速度と位置を修正して
            このタイムステップの値を確定する。
            (位置の修正量は a × Δt^2)
            使い終わった加速度は全粒子について 0 に戻す。
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】Velocity[], Position[], Acceleration[]
  【呼び出し元】mainLoopOfSimulation()
=====================================================================*/
void moveParticleUsingPressureGradient( void ){
  int i;   /* 粒子番号 */

  for(i=0;i<NumberOfParticles;i++){
    if(ParticleType[i] == FLUID){
      Velocity[i*3  ] +=Acceleration[i*3  ]*DT;
      Velocity[i*3+1] +=Acceleration[i*3+1]*DT;
      Velocity[i*3+2] +=Acceleration[i*3+2]*DT;

      Position[i*3  ] +=Acceleration[i*3  ]*DT*DT;
      Position[i*3+1] +=Acceleration[i*3+1]*DT*DT;
      Position[i*3+2] +=Acceleration[i*3+2]*DT*DT;
    }
    Acceleration[i*3  ]=0.0;
    Acceleration[i*3+1]=0.0;
    Acceleration[i*3+2]=0.0;
  }
}


/*=====================================================================
  【関数名】writeData_inProfFormat
  【機能】  現在の状態を output_%04d.prof (テキスト形式) に書き出し、
            ファイル番号を 1 進める。
            1 行目: 時刻, 2 行目: 粒子数,
            3 行目以降: 種類 x y z vx vy vz 圧力 粒子数密度 (1 粒子 1 行)
  【引数】  なし (void)
  【戻り値】なし (void)
  【書き換えるグローバル変数】FileNumber
  【呼び出し元】mainLoopOfSimulation()
  【注意】  writeData_inVtuFormat() の後に呼ぶこと (番号を進めるため)
=====================================================================*/
void writeData_inProfFormat( void ){
  int i;               /* 粒子番号 */
  FILE *fp;            /* 出力ファイルのポインタ */
  char fileName[256];  /* 出力ファイル名 */

  sprintf(fileName, "output_%04d.prof",FileNumber);
  fp = fopen(fileName, "w");
  fprintf(fp,"%lf\n",Time);
  fprintf(fp,"%d\n",NumberOfParticles);
  for(i=0;i<NumberOfParticles;i++) {
    fprintf(fp,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n"
	    ,ParticleType[i], Position[i*3], Position[i*3+1], Position[i*3+2]
	    ,Velocity[i*3], Velocity[i*3+1], Velocity[i*3+2], Pressure[i], NumberDensity[i]);
  }
  fclose(fp);
  FileNumber++;
}


/*=====================================================================
  【関数名】writeData_inVtuFormat
  【機能】  現在の状態を particle_%04d.vtu (VTK の UnstructuredGrid,
            XML 形式) に書き出す。各粒子を 1 点のセル (VTK_VERTEX) とし、
            粒子の種類・速さ・圧力を点データとして出力する。
            ParaView などでそのまま開ける。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】mainLoopOfSimulation()
  【注意】  ファイル番号は進めない (writeData_inProfFormat() が進める)
=====================================================================*/
void writeData_inVtuFormat( void ){
  int i;                           /* 粒子番号 */
  double absoluteValueOfVelocity;  /* 速度の大きさ (速さ) [m/s] */
  FILE *fp;                        /* 出力ファイルのポインタ */
  char fileName[1024];             /* 出力ファイル名 */

  sprintf(fileName, "particle_%04d.vtu", FileNumber);
  fp=fopen(fileName,"w");
  fprintf(fp,"<?xml version='1.0' encoding='UTF-8'?>\n");
  fprintf(fp,"<VTKFile xmlns='VTK' byte_order='LittleEndian' version='0.1' type='UnstructuredGrid'>\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfCells='%d' NumberOfPoints='%d'>\n",NumberOfParticles,NumberOfParticles);
  /* 粒子の位置 */
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray NumberOfComponents='3' type='Float32' Name='Position' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%lf %lf %lf\n",Position[i*3],Position[i*3+1],Position[i*3+2]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  /* 粒子ごとのデータ (種類・速さ・圧力) */
  fprintf(fp,"<PointData>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Int32' Name='ParticleType' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",ParticleType[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='Velocity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    absoluteValueOfVelocity=
      sqrt( Velocity[i*3]*Velocity[i*3] + Velocity[i*3+1]*Velocity[i*3+1] + Velocity[i*3+2]*Velocity[i*3+2] );
    fprintf(fp,"%f\n",(float)absoluteValueOfVelocity);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray NumberOfComponents='1' type='Float32' Name='pressure' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%f\n",(float)Pressure[i]);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");
  /* セル情報 (粒子 1 個 = 点 1 個のセル) */
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type='Int32' Name='connectivity' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='Int32' Name='offsets' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"%d\n",i+1);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type='UInt8' Name='types' format='ascii'>\n");
  for(i=0;i<NumberOfParticles;i++){
    fprintf(fp,"1\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");
  fprintf(fp,"</Piece>\n");
  fprintf(fp,"</UnstructuredGrid>\n");
  fprintf(fp,"</VTKFile>\n");
  fclose(fp);
}
