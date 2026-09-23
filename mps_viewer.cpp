/*=====================================================================
  mps_viewer.cpp

   MPS 法シミュレーション (mps.c) の計算結果を
   DirectX 11 でアニメーション表示するビューア

   (c) このファイルは mps.c (Kazuya SHIBATA, Kohei MUROTANI and
       Seiichi KOSHIZUKA, 2014) の結果を可視化するために作成した。
=======================================================================*/
/*=====================================================================
 【このプログラムの概要】

   mps.c が出力する output_%04d.prof を全部読み込んでメモリに置き、
   DirectX 11 で 1 粒子 = 1 個の円 (ビルボード) として描画し、
   パラパラ漫画の要領でアニメーション再生する。

 【全体の流れ】

     main()
       ├ ParseCommandLine()  : コマンドライン引数の解釈
       ├ LoadAllFrames()     : output_*.prof を全フレーム読み込む
       ├ AnalyzeFrames()     : 粒子間距離・値の範囲・境界箱を調べる
       ├ CreateAppWindow()   : Win32 ウィンドウを作る
       ├ InitD3D()           : DirectX 11 の初期化
       │    ├ D3D11CreateDeviceAndSwapChain() : デバイスとスワップチェーン
       │    ├ CreateRenderTargets()           : 描画先とZバッファ
       │    ├ CompileShaders()                : HLSL を実行時コンパイル
       │    └ CreateBuffers()                 : 頂点/インスタンス/定数バッファ
       ├ MainLoop()          : メッセージ処理 + 毎フレーム描画
       │    ├ UpdateAnimation() : 経過時間からフレーム番号を進める
       │    ├ BuildInstances()  : 粒子ごとの位置と色を作る
       │    └ RenderFrame()     : GPU に描画させて画面に出す
       └ Cleanup()           : 後片付け

 【描画方式 (インスタンシング)】

   粒子は 1 個ずつ「カメラの方を向いた正方形 (ビルボード)」として
   描く。正方形の 4 頂点は 1 個だけ用意し、粒子ごとに異なるのは
   「中心位置」と「色」だけなので、DrawInstanced() で
   粒子数分まとめて描画する。ピクセルシェーダで正方形の角を
   切り落として円形にし、簡易的な陰影を付けて球のように見せる。

   ・頂点バッファ       : 正方形の 4 隅 (全粒子で共有)
   ・インスタンスバッファ: 粒子ごとの [中心座標 xyz][色 rgba]
   ・定数バッファ       : ビュー射影行列、ビルボードの右/上ベクトル

 【必要な環境】
   ・Windows 10 / 11
   ・Visual Studio 2022 (Windows SDK 10 に DirectX 11 が含まれる)
   ・レガシーの DirectX SDK (June 2010) は不要

 【ビルド方法】
   ● 開発者コマンドプロンプトで
         cl /W3 /O2 /EHsc mps_viewer.cpp /Fe:mps_viewer.exe
     (必要なライブラリは #pragma comment で自動リンクされる)

   ● Visual Studio の IDE で
         「空のプロジェクト」(C++) を作り、本ファイルを追加して
         Release / x64 でビルドする。

 【実行方法】
         mps_viewer.exe [結果フォルダ] [オプション]

     結果フォルダ : output_0000.prof などが置いてあるフォルダ。
                    省略するとカレントディレクトリを見る。

     オプション
       -r <半径>          粒子の描画半径 [m] を手動指定
                          (省略時は粒子間距離から自動推定)
       -c <0-3>           起動時の配色モード
                          0:圧力 1:速度 2:粒子種類 3:粒子数密度
       -wall / -nowall    壁粒子を表示する / しない
                          (省略時は 2 次元なら表示、3 次元なら非表示。
                           3 次元では壁が閉じた箱になっており、表示すると
                           外側の壁しか見えず中の流体が確認できないため)
       -yaw   <角度>      起動時のカメラ水平回転角 [度]（3 次元計算用）
       -pitch <角度>      起動時のカメラ仰角 [度]（3 次元計算用）
                          例: -yaw 35 -pitch 20 で斜め上から見下ろす
       -shot <番号> <ファイル名.bmp>
                          画面を出さずに指定フレームを BMP に保存して終了
                          (動作確認やレポート用の静止画作成に使う)

 【操作方法】(ウィンドウをアクティブにして操作する)
     Space        再生 / 一時停止
     → / ←       1 フレーム進む / 戻る
     Home         先頭フレームへ
     L            ループ再生の ON/OFF
     + / -        再生速度を上げる / 下げる
     0            再生速度を等倍に戻す
     1 / 2 / 3 / 4  配色モード (圧力 / 速度 / 粒子種類 / 粒子数密度)
     [ / ]        カラースケールの上限を下げる / 上げる
     W            壁粒子の表示 ON/OFF
     D            ダミー壁粒子の表示 ON/OFF
     , / .        粒子の描画半径を小さく / 大きく
     F            表示範囲を全体が入るように戻す (フィット)
     F12          スクリーンショットを shot_%03d.bmp に保存
     Esc          終了

     マウス左ドラッグ : 回転 (3 次元計算の結果を見るとき)
     マウス右ドラッグ : 平行移動
     ホイール         : 拡大 / 縮小

   ※画面に文字を描く機能は持たせず、現在のフレーム番号・時刻・
     配色モードなどの情報はウィンドウのタイトルバーに表示している。
     文字描画のために DirectWrite を組み込むと本質でないコードが
     増えるため、あえて簡素にした。
=======================================================================*/

#define _CRT_SECURE_NO_WARNINGS   /* fopen/sprintf の警告 C4996 を抑止 */
#define WIN32_LEAN_AND_MEAN       /* windows.h の読み込みを軽くする    */
#define NOMINMAX                  /* min/max マクロを無効化 (std と衝突) */

#include <windows.h>
#include <d3d11.h>                /* DirectX 11 本体                   */
#include <d3dcompiler.h>          /* HLSL の実行時コンパイル           */
#include <DirectXMath.h>          /* 行列・ベクトル演算                */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <string>
#include <algorithm>   /* std::nth_element 用 */

/* 必要なライブラリを自動でリンクする (VC++ 専用の書き方) */
#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3dcompiler.lib")
#pragma comment(lib, "user32.lib")      /* ウィンドウ関連の API */

using namespace DirectX;

/*---------------------------------------------------------------------
  マウス座標を LPARAM から取り出すマクロ。

  windowsx.h の GET_X_LPARAM と同じもの。
  単純に LOWORD() を使うと、マウスをウィンドウの外へドラッグした
  ときの負の座標が巨大な正の値になってしまうため、
  いったん short (符号付き 16bit) にキャストしてから int にする。
---------------------------------------------------------------------*/
#define GET_X_LPARAM_COMPAT(lp) ((int)(short)LOWORD(lp))
#define GET_Y_LPARAM_COMPAT(lp) ((int)(short)HIWORD(lp))

/*---------------------------------------------------------------------
  COM オブジェクトを安全に解放するためのマクロ。
  DirectX のオブジェクトは Release() で参照カウントを減らして解放する。
---------------------------------------------------------------------*/
#define SAFE_RELEASE(p) do { if(p){ (p)->Release(); (p) = nullptr; } } while(0)

/*---------------------------------------------------------------------
  【粒子種類】mps.c の #define と同じ値を使う
---------------------------------------------------------------------*/
enum ParticleTypeEnum {
	PT_GHOST      = -1,   /* 計算対象外           */
	PT_FLUID      =  0,   /* 流体粒子             */
	PT_WALL       =  2,   /* 壁粒子               */
	PT_DUMMY_WALL =  3    /* ダミー壁粒子         */
};

/*---------------------------------------------------------------------
  【配色モード】どの物理量で粒子に色を付けるか
---------------------------------------------------------------------*/
enum ColorMode {
	CM_PRESSURE = 0,      /* 圧力         [Pa]   */
	CM_VELOCITY = 1,      /* 速度の大きさ [m/s]  */
	CM_TYPE     = 2,      /* 粒子の種類          */
	CM_DENSITY  = 3,      /* 粒子数密度   [-]    */
	CM_COUNT    = 4
};
static const char* kColorModeName[CM_COUNT] = { "Pressure", "Velocity", "Type", "NumDensity" };

/*=====================================================================
  【構造体】読み込んだシミュレーション結果を保持する
=====================================================================*/

/* 粒子 1 個分のデータ (.prof の 1 行に対応する) */
struct Particle {
	float px, py, pz;     /* 位置         [m]    */
	float vx, vy, vz;     /* 速度         [m/s]  */
	float pressure;       /* 圧力         [Pa]   */
	float numberDensity;  /* 粒子数密度   [-]    */
	float speed;          /* 速度の大きさ [m/s] (読み込み時に計算しておく) */
	int   type;           /* 粒子の種類 (ParticleTypeEnum)                 */
};

/* 1 時刻分のデータ (.prof ファイル 1 個に対応する) */
struct Frame {
	float time;                      /* そのフレームの計算時刻 [s] */
	std::vector<Particle> particles; /* 全粒子のデータ             */
};

/* GPU へ送る「粒子 1 個分」のインスタンスデータ。
   HLSL 側の INSTPOS / INSTCOL と並び順を一致させること。 */
struct InstanceData {
	float x, y, z;        /* 粒子の中心座標 (ワールド座標)  */
	float r, g, b, a;     /* 描画色 (RGBA, 0.0～1.0)        */
};

/* GPU へ送る定数バッファ。
   HLSL の cbuffer と並び・サイズを一致させること。
   DirectX の定数バッファは 16 バイト単位でないといけないため、
   float3 ではなく float4 を使って詰め物 (padding) を入れている。 */
struct ConstantBufferData {
	XMFLOAT4X4 viewProj;  /* ビュー行列 × 射影行列 (転置済み)      */
	XMFLOAT4   right;     /* ビルボードの右方向ベクトル × 粒子半径 */
	XMFLOAT4   up;        /* ビルボードの上方向ベクトル × 粒子半径 */
};

/*=====================================================================
  【グローバル変数】
  mps.c と同様、状態はグローバル変数で持ち、関数間で共有する。
=====================================================================*/

/* --- シミュレーション結果 --- */
static std::vector<Frame> g_frames;        /* 全フレームのデータ        */
static int    g_maxParticles   = 0;        /* 全フレーム中の最大粒子数  */
static float  g_frameDeltaTime = 0.02f;    /* フレーム間の計算時間差 [s]*/
static float  g_particleSpacing = 0.025f;  /* 推定した粒子間距離 l0 [m] */
static bool   g_isThreeDimensional = false;
                        /* 読み込んだデータが 3 次元計算かどうか。
                           z 方向に広がりがあれば 3 次元と判定する。
                           3 次元では壁が閉じた箱になり、外から見ると
                           壁しか見えないため、既定で壁を非表示にする */

/* 物理量の最小値・最大値 (カラースケールの基準に使う) */
static float  g_pressureMin = 0.0f, g_pressureMax = 1.0f;
static float  g_speedMin    = 0.0f, g_speedMax    = 1.0f;
static float  g_densityMin  = 0.0f, g_densityMax  = 1.0f;

/* 全粒子を包む境界箱 (カメラの初期位置を決めるのに使う) */
static XMFLOAT3 g_boundsMin = { 0, 0, 0 };
static XMFLOAT3 g_boundsMax = { 1, 1, 0 };

/* --- 再生状態 --- */
static float  g_currentFrame   = 0.0f;   /* 現在のフレーム番号 (小数)  */
static bool   g_isPlaying      = true;   /* 再生中かどうか             */
static bool   g_isLooping      = true;   /* 末尾で先頭へ戻るか         */
static float  g_playbackSpeed  = 1.0f;   /* 再生速度の倍率             */
static int    g_colorMode      = CM_PRESSURE;  /* 配色モード           */
static float  g_colorScale     = 1.0f;   /* カラースケール上限の倍率   */
static bool   g_showWall       = true;   /* 壁粒子を表示するか         */
static bool   g_showDummyWall  = false;  /* ダミー壁粒子を表示するか   */
static float  g_particleRadius = 0.0f;   /* 粒子の描画半径 [m]
                                            (0 なら自動推定)           */
static int    g_screenshotCount = 0;     /* F12 で保存した枚数         */

/* --- カメラ (正射影 + 軌道 (オービット) 操作) --- */
static float    g_cameraYaw      = 0.0f;  /* 水平回転角 [rad]          */
static float    g_cameraPitch    = 0.0f;  /* 仰角       [rad]          */
static XMFLOAT3 g_cameraTarget   = { 0, 0, 0 }; /* 注視点             */
static float    g_cameraDistance = 10.0f; /* 注視点からの距離 [m]      */
static float    g_orthoHeight    = 1.0f;  /* 画面に映る高さ [m]
                                             小さいほど拡大される      */

/* --- マウス操作の状態 --- */
static bool  g_isDraggingLeft  = false;   /* 左ボタンで回転中          */
static bool  g_isDraggingRight = false;   /* 右ボタンで平行移動中      */
static POINT g_lastMousePos    = { 0, 0 };

/* --- Win32 / DirectX のオブジェクト --- */
static HWND                     g_hWnd            = nullptr;
static ID3D11Device*            g_device          = nullptr; /* リソース生成役 */
static ID3D11DeviceContext*     g_context         = nullptr; /* 描画命令発行役 */
static IDXGISwapChain*          g_swapChain       = nullptr; /* 画面への出力   */
static ID3D11RenderTargetView*  g_renderTargetView= nullptr; /* 描画先 (色)    */
static ID3D11Texture2D*         g_depthStencil    = nullptr; /* Z バッファ実体 */
static ID3D11DepthStencilView*  g_depthStencilView= nullptr; /* Z バッファ参照 */
static ID3D11VertexShader*      g_vertexShader    = nullptr;
static ID3D11PixelShader*       g_pixelShader     = nullptr;
static ID3D11InputLayout*       g_inputLayout     = nullptr; /* 頂点の並び定義 */
static ID3D11Buffer*            g_quadVertexBuffer= nullptr; /* 正方形 4 頂点  */
static ID3D11Buffer*            g_instanceBuffer  = nullptr; /* 粒子ごとの情報 */
static ID3D11Buffer*            g_constantBuffer  = nullptr; /* 行列など       */
static ID3D11RasterizerState*   g_rasterizerState = nullptr;
static ID3D11DepthStencilState* g_depthState      = nullptr;

static int g_backBufferWidth  = 1280;     /* 描画先の幅  [pixel]       */
static int g_backBufferHeight = 720;      /* 描画先の高さ[pixel]       */

/* インスタンスデータを組み立てる作業用配列 (毎フレーム使い回す) */
static std::vector<InstanceData> g_instances;

/*=====================================================================
  【HLSL シェーダのソース】

  実行時に D3DCompile() でコンパイルするため、文字列として埋め込む。
  こうしておくと .hlsl ファイルや fxc.exe によるビルド手順が不要になり、
  exe 1 個だけで動くようになる。

  ●頂点シェーダ VSMain
     粒子の中心座標に、カメラの右方向 × corner.x と
     上方向 × corner.y を足して正方形の 4 隅を作る。
     こうすると正方形が常にカメラの方を向く (ビルボード)。

  ●ピクセルシェーダ PSMain
     正方形の中心からの距離が 1 を超えるピクセルを clip() で捨てて
     円形にする。さらに中心ほど明るくして球のように見せる。
=====================================================================*/
static const char* kShaderSource = R"HLSL(
cbuffer ConstantBuffer : register(b0)
{
    float4x4 gViewProj;   // ビュー行列 × 射影行列
    float4   gRight;      // カメラの右方向 × 粒子半径
    float4   gUp;         // カメラの上方向 × 粒子半径
};

struct VSInput
{
    float2 corner   : POSITION;   // 頂点ごと  : 正方形の隅 (-1～+1)
    float3 instPos  : INSTPOS;    // 粒子ごと  : 中心座標
    float4 instCol  : INSTCOL;    // 粒子ごと  : 色
};

struct VSOutput
{
    float4 position : SV_Position; // 画面上の座標 (必須)
    float4 color    : COLOR;       // 色
    float2 uv       : TEXCOORD0;   // 正方形内の位置 (-1～+1)
};

VSOutput VSMain(VSInput input)
{
    // 粒子の中心から、カメラの右・上方向へ広げて正方形の頂点を作る
    float3 worldPos = input.instPos
                    + gRight.xyz * input.corner.x
                    + gUp.xyz    * input.corner.y;

    VSOutput output;
    output.position = mul(float4(worldPos, 1.0f), gViewProj);
    output.color    = input.instCol;
    output.uv       = input.corner;
    return output;
}

float4 PSMain(VSOutput input) : SV_Target
{
    // 中心からの距離の2乗。1 を超える = 正方形の角なので描かない
    float r2 = dot(input.uv, input.uv);
    clip(1.0f - r2);

    // 中心ほど明るくして球のような陰影を付ける
    float shade = sqrt(saturate(1.0f - r2));
    float3 rgb  = input.color.rgb * (0.40f + 0.60f * shade);

    // 輪郭をわずかに暗くして粒子どうしの境界を見やすくする
    rgb *= lerp(1.0f, 0.75f, saturate((r2 - 0.55f) / 0.45f));

    return float4(rgb, 1.0f);
}
)HLSL";

/*=====================================================================
  【関数名】Clampf
  【機能】  値 v を [lo, hi] の範囲に収める (クランプする)。
  【引数】  float v  : 対象の値
            float lo : 下限
            float hi : 上限
  【戻り値】float : 範囲内に収めた値
=====================================================================*/
static float Clampf(float v, float lo, float hi)
{
	if (v < lo) return lo;
	if (v > hi) return hi;
	return v;
}

/*=====================================================================
  【関数名】PercentileOf
  【機能】  値の配列から「下から ratio の割合の位置にある値」
            (パーセンタイル) を求める。

            なぜ最大値ではなくパーセンタイルを使うのか:
            MPS 法の圧力は、粒子が瞬間的にぶつかった箇所で
            ごく一部だけ極端に大きな値 (スパイク) になることがある。
            その最大値をカラースケールの上限にすると、
            他の粒子はすべて下端の色 (青) に潰れて何も見えなくなる。
            そこで上位 1% を外れ値として切り捨て、99 パーセンタイルを
            上限にすることで、全体の圧力分布が見えるようにする。

  【引数】  std::vector<float>& values : 対象の値 (並べ替えられる)
            float ratio               : 0.0〜1.0 の割合 (0.99 なら99%)
  【戻り値】float : パーセンタイル値。配列が空なら 0
  【呼び出し元】AnalyzeFrames()
=====================================================================*/
static float PercentileOf(std::vector<float>& values, float ratio)
{
	if (values.empty()) return 0.0f;

	size_t index = (size_t)(ratio * (float)(values.size() - 1));
	if (index >= values.size()) index = values.size() - 1;

	/* 全体を並べ替える必要はなく、index 番目だけ確定すればよいので
	   std::sort より速い std::nth_element を使う */
	std::nth_element(values.begin(), values.begin() + index, values.end());
	return values[index];
}

/*=====================================================================
  【関数名】JetColorMap
  【機能】  0.0～1.0 の値を「青 → 水色 → 緑 → 黄 → 赤」の色に変換する。
            数値可視化で広く使われる jet カラーマップの簡易版。
            低い値ほど青く、高い値ほど赤くなる。
  【引数】  float t   : 正規化した値 (0.0～1.0)
            float* rgb: 結果を受け取る配列 (要素 3 個)
  【戻り値】なし (void)。結果は引数 rgb に書き込む。
  【呼び出し元】BuildInstances()
=====================================================================*/
static void JetColorMap(float t, float* rgb)
{
	t = Clampf(t, 0.0f, 1.0f);

	/* 0.00 青(0,0,1) → 0.25 水(0,1,1) → 0.50 緑(0,1,0)
	   → 0.75 黄(1,1,0) → 1.00 赤(1,0,0)                 */
	if (t < 0.25f) {
		float s = t / 0.25f;
		rgb[0] = 0.0f;  rgb[1] = s;     rgb[2] = 1.0f;
	} else if (t < 0.50f) {
		float s = (t - 0.25f) / 0.25f;
		rgb[0] = 0.0f;  rgb[1] = 1.0f;  rgb[2] = 1.0f - s;
	} else if (t < 0.75f) {
		float s = (t - 0.50f) / 0.25f;
		rgb[0] = s;     rgb[1] = 1.0f;  rgb[2] = 0.0f;
	} else {
		float s = (t - 0.75f) / 0.25f;
		rgb[0] = 1.0f;  rgb[1] = 1.0f - s; rgb[2] = 0.0f;
	}
}

/*=====================================================================
  【関数名】LoadOneFrame
  【機能】  .prof ファイルを 1 個読み込んで Frame 構造体に格納する。

            .prof の書式 (mps.c の writeData_inProfFormat() が出力)
              1 行目   : 時刻
              2 行目   : 粒子数
              3 行目以降: 種類 x y z vx vy vz 圧力 粒子数密度

  【引数】  const char* path : 読み込むファイルのパス
            Frame& frame     : 結果を格納する構造体 (参照渡し)
  【戻り値】bool : 読み込めたら true、ファイルが無ければ false
  【呼び出し元】LoadAllFrames()
=====================================================================*/
static bool LoadOneFrame(const char* path, Frame& frame)
{
	FILE* fp = fopen(path, "r");
	if (fp == nullptr) return false;   /* ここで終端を検出する */

	int numParticles = 0;
	if (fscanf(fp, "%f", &frame.time) != 1 ||
	    fscanf(fp, "%d", &numParticles) != 1 || numParticles <= 0) {
		fclose(fp);
		return false;
	}

	frame.particles.clear();
	frame.particles.reserve(numParticles);

	for (int i = 0; i < numParticles; i++) {
		Particle p;
		int n = fscanf(fp, "%d %f %f %f %f %f %f %f %f",
		               &p.type,
		               &p.px, &p.py, &p.pz,
		               &p.vx, &p.vy, &p.vz,
		               &p.pressure, &p.numberDensity);
		if (n != 9) break;   /* 途中で壊れていたらそこまでを使う */

		/* 速度の大きさは毎フレーム使うので、ここで計算しておく */
		p.speed = sqrtf(p.vx * p.vx + p.vy * p.vy + p.vz * p.vz);
		frame.particles.push_back(p);
	}

	fclose(fp);
	return !frame.particles.empty();
}

/*=====================================================================
  【関数名】LoadAllFrames
  【機能】  指定フォルダの output_0000.prof から順に、ファイルが
            見つからなくなるまで連続して読み込む。
  【引数】  const char* directory : .prof が置いてあるフォルダ
  【戻り値】bool : 1 フレーム以上読み込めたら true
  【副作用】g_frames に全フレームを格納する
  【呼び出し元】main()
=====================================================================*/
static bool LoadAllFrames(const char* directory)
{
	g_frames.clear();

	char path[1024];
	for (int index = 0; ; index++) {
		sprintf(path, "%s\\output_%04d.prof", directory, index);

		Frame frame;
		if (!LoadOneFrame(path, frame)) break;   /* 見つからない = 終端 */

		g_frames.push_back(std::move(frame));

		/* 大量に読み込むときのために進捗を出す */
		if ((index % 50) == 0) {
			printf("  loading ... output_%04d.prof (%d particles)\n",
			       index, (int)g_frames.back().particles.size());
		}
	}

	return !g_frames.empty();
}

/*=====================================================================
  【関数名】AnalyzeFrames
  【機能】  読み込んだ全フレームを走査して、描画に必要な情報を調べる。

              ・最大粒子数         → GPU バッファの確保サイズ
              ・境界箱 (bounding box) → カメラの初期位置
              ・圧力/速度/粒子数密度の最小・最大値 → カラースケール
              ・粒子間距離 l0      → 粒子の描画半径
              ・フレーム間の時間差 → 再生速度

            粒子間距離は、先頭フレームの粒子から最大 400 個を抜き出し、
            総当たりで最も近い 2 粒子の距離を求めて推定する。
            MPS 法では粒子が格子状に並ぶので、この値がほぼ
            PARTICLE_DISTANCE (l0) に一致する。

  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】g_maxParticles, g_bounds*, g_pressure*, g_speed*,
            g_density*, g_particleSpacing, g_frameDeltaTime を設定
  【呼び出し元】main()
=====================================================================*/
static void AnalyzeFrames(void)
{
	g_maxParticles = 0;

	/* カラースケールをパーセンタイルで決めるため、
	   流体粒子の値をいったん全部集めてから統計をとる */
	std::vector<float> pressureValues;
	std::vector<float> speedValues;
	std::vector<float> densityValues;

	for (size_t f = 0; f < g_frames.size(); f++) {
		const Frame& frame = g_frames[f];
		if ((int)frame.particles.size() > g_maxParticles) {
			g_maxParticles = (int)frame.particles.size();
		}
		for (size_t i = 0; i < frame.particles.size(); i++) {
			const Particle& p = frame.particles[i];
			/* カラースケールは流体粒子の値で決める
			   (壁やダミー壁を含めると範囲が偏るため) */
			if (p.type == PT_FLUID) {
				pressureValues.push_back(p.pressure);
				speedValues.push_back(p.speed);
				densityValues.push_back(p.numberDensity);
			}
		}
	}

	/* --- 境界箱 (表示範囲) は「先頭フレーム」だけから決める ---
	   全フレームの最大範囲で合わせると、飛び散った数個のしぶきに
	   引っ張られて水槽がごく小さくしか表示されなくなってしまう。
	   初期配置 = 水槽そのものなので、これを基準にするのが最も自然。
	   ダミー壁は水槽の外側に広がっていて邪魔なので除外する。       */
	{
		float xMin =  1e30f, yMin =  1e30f, zMin =  1e30f;
		float xMax = -1e30f, yMax = -1e30f, zMax = -1e30f;

		const std::vector<Particle>& ps = g_frames[0].particles;
		for (size_t i = 0; i < ps.size(); i++) {
			const Particle& p = ps[i];
			if (p.type == PT_DUMMY_WALL || p.type == PT_GHOST) continue;
			if (p.px < xMin) xMin = p.px;   if (p.px > xMax) xMax = p.px;
			if (p.py < yMin) yMin = p.py;   if (p.py > yMax) yMax = p.py;
			if (p.pz < zMin) zMin = p.pz;   if (p.pz > zMax) zMax = p.pz;
		}
		if (xMin > xMax) { xMin = yMin = zMin = 0.0f; xMax = yMax = zMax = 1.0f; }

		g_boundsMin = XMFLOAT3(xMin, yMin, zMin);
		g_boundsMax = XMFLOAT3(xMax, yMax, zMax);
	}

	/* --- カラースケールの範囲 ---
	   上限は最大値ではなく 99 パーセンタイルを使う。
	   圧力は衝突の瞬間に一部の粒子だけ極端な値になるため、
	   最大値を上限にすると全体が青一色に潰れてしまう。         */
	{
		float pMin = pressureValues.empty() ? 0.0f : *std::min_element(pressureValues.begin(), pressureValues.end());
		float sMin = speedValues.empty()    ? 0.0f : *std::min_element(speedValues.begin(),    speedValues.end());
		float dMin = densityValues.empty()  ? 0.0f : *std::min_element(densityValues.begin(),  densityValues.end());

		float pMax = PercentileOf(pressureValues, 0.99f);
		float sMax = PercentileOf(speedValues,    0.99f);
		float dMax = densityValues.empty() ? 1.0f
		           : *std::max_element(densityValues.begin(), densityValues.end());
		/* 粒子数密度は外れ値が出にくいので最大値をそのまま使う */

		g_pressureMin = pMin;  g_pressureMax = (pMax > pMin) ? pMax : pMin + 1.0f;
		g_speedMin    = sMin;  g_speedMax    = (sMax > sMin) ? sMax : sMin + 1.0f;
		g_densityMin  = dMin;  g_densityMax  = (dMax > dMin) ? dMax : dMin + 1.0f;
	}

	/* --- 粒子間距離 l0 の推定 (最近接粒子間距離を総当たりで探す) --- */
	{
		const std::vector<Particle>& ps = g_frames[0].particles;
		int sampleCount = (int)ps.size();
		if (sampleCount > 400) sampleCount = 400;

		float minDist2 = 1e30f;
		for (int i = 0; i < sampleCount; i++) {
			for (int j = i + 1; j < sampleCount; j++) {
				float dx = ps[i].px - ps[j].px;
				float dy = ps[i].py - ps[j].py;
				float dz = ps[i].pz - ps[j].pz;
				float d2 = dx * dx + dy * dy + dz * dz;
				if (d2 > 1e-18f && d2 < minDist2) minDist2 = d2;
			}
		}
		if (minDist2 < 1e29f) g_particleSpacing = sqrtf(minDist2);
	}

	/* --- フレーム間の計算時間差 --- */
	if (g_frames.size() >= 2) {
		float dt = g_frames[1].time - g_frames[0].time;
		if (dt > 1e-9f) g_frameDeltaTime = dt;
	}

	/* 描画半径が未指定なら粒子間距離の半分にする
	   (粒子どうしがちょうど接するように見える) */
	if (g_particleRadius <= 0.0f) {
		g_particleRadius = 0.5f * g_particleSpacing;
	}

	/* --- 2 次元計算か 3 次元計算かを判定する ---
	   2 次元計算では全粒子の z 座標が 0 なので z 方向の広がりは 0 になる。
	   粒子間距離の半分より広がっていれば 3 次元とみなす。           */
	g_isThreeDimensional =
		((g_boundsMax.z - g_boundsMin.z) > 0.5f * g_particleSpacing);
}

/*=====================================================================
  【関数名】FitCameraToBounds
  【機能】  境界箱の全体が画面に収まるようにカメラを設定し直す。
  【引数】  なし (void)
  【戻り値】なし (void)
  【副作用】g_cameraTarget, g_orthoHeight, g_cameraDistance を設定
  【呼び出し元】main() (起動時), WndProc() (F キー)
=====================================================================*/
static void FitCameraToBounds(void)
{
	float cx = 0.5f * (g_boundsMin.x + g_boundsMax.x);
	float cy = 0.5f * (g_boundsMin.y + g_boundsMax.y);
	float cz = 0.5f * (g_boundsMin.z + g_boundsMax.z);
	g_cameraTarget = XMFLOAT3(cx, cy, cz);

	float sx = g_boundsMax.x - g_boundsMin.x;
	float sy = g_boundsMax.y - g_boundsMin.y;
	float sz = g_boundsMax.z - g_boundsMin.z;

	/* 画面の縦横比を考慮して、横にも縦にもはみ出さない表示高さを求める。
	   正射影では「画面に映る高さ」だけを指定し、幅は
	   高さ × 縦横比 で決まる。したがって横方向に必要な大きさは
	   縦横比で割って高さに換算してから比較する。
	   これをしないと横長の計算領域が画面の中央に小さく表示される。 */
	float aspect = (g_backBufferHeight > 0)
	             ? (float)g_backBufferWidth / (float)g_backBufferHeight : 1.0f;

	/* 回転させて奥行き方向を向いたときも収まるよう、
	   横方向は x と z の大きいほうを使う */
	float horizontal = (sx > sz) ? sx : sz;
	float vertical   = sy;

	float requiredHeight = vertical;
	if (horizontal / aspect > requiredHeight) requiredHeight = horizontal / aspect;
	if (requiredHeight < 1e-6f) requiredHeight = 1.0f;

	g_orthoHeight = requiredHeight * 1.10f;   /* 少し余白を付ける */

	/* 正射影では視点の距離は見た目の大きさに影響しない。
	   Z バッファの範囲に全粒子が入るよう十分遠くに置く */
	float diagonal = sqrtf(sx * sx + sy * sy + sz * sz);
	g_cameraDistance = (diagonal > 1e-6f) ? diagonal * 8.0f : 10.0f;
}

/*=====================================================================
  【関数名】BuildInstances
  【機能】  現在のフレームの粒子から、GPU へ送るインスタンスデータ
            (中心座標 + 色) を作る。
            表示が OFF の粒子種類はここで除外するので、
            GPU へ送る個数は毎フレーム変わる。
  【引数】  int frameIndex : 描画するフレーム番号
  【戻り値】int : 実際に描画する粒子の個数
  【副作用】g_instances に結果を詰める
  【呼び出し元】RenderFrame()
=====================================================================*/
static int BuildInstances(int frameIndex)
{
	g_instances.clear();

	if (frameIndex < 0 || frameIndex >= (int)g_frames.size()) return 0;
	const std::vector<Particle>& ps = g_frames[frameIndex].particles;

	/* 配色モードごとの値の範囲を決める。
	   g_colorScale は [ ] キーで変えられる上限の倍率 */
	float valueMin = 0.0f, valueMax = 1.0f;
	if (g_colorMode == CM_PRESSURE) {
		valueMin = g_pressureMin;
		valueMax = g_pressureMin + (g_pressureMax - g_pressureMin) * g_colorScale;
	} else if (g_colorMode == CM_VELOCITY) {
		valueMin = g_speedMin;
		valueMax = g_speedMin + (g_speedMax - g_speedMin) * g_colorScale;
	} else if (g_colorMode == CM_DENSITY) {
		valueMin = g_densityMin;
		valueMax = g_densityMin + (g_densityMax - g_densityMin) * g_colorScale;
	}
	float valueRange = valueMax - valueMin;
	if (valueRange < 1e-12f) valueRange = 1e-12f;

	for (size_t i = 0; i < ps.size(); i++) {
		const Particle& p = ps[i];

		/* --- 表示するかどうかの判定 --- */
		if (p.type == PT_GHOST) continue;
		if (p.type == PT_WALL       && !g_showWall)      continue;
		if (p.type == PT_DUMMY_WALL && !g_showDummyWall) continue;

		InstanceData inst;
		inst.x = p.px;  inst.y = p.py;  inst.z = p.pz;
		inst.a = 1.0f;

		if (g_colorMode == CM_TYPE) {
			/* 粒子の種類で色分けする */
			if (p.type == PT_FLUID) {
				inst.r = 0.25f; inst.g = 0.55f; inst.b = 1.00f;  /* 青  */
			} else if (p.type == PT_WALL) {
				inst.r = 0.70f; inst.g = 0.70f; inst.b = 0.72f;  /* 灰  */
			} else {
				inst.r = 0.35f; inst.g = 0.33f; inst.b = 0.38f;  /* 暗灰 */
			}
		} else if (p.type == PT_DUMMY_WALL) {
			/* ダミー壁は圧力を解いていない (常に 0) ので、
			   物理量で色を付けても意味がない。常に暗い灰色にする */
			inst.r = 0.30f; inst.g = 0.29f; inst.b = 0.33f;
		} else {
			/* 物理量を 0～1 に正規化してカラーマップに通す */
			float value = 0.0f;
			if      (g_colorMode == CM_PRESSURE) value = p.pressure;
			else if (g_colorMode == CM_VELOCITY) value = p.speed;
			else                                 value = p.numberDensity;

			float t = (value - valueMin) / valueRange;
			float rgb[3];
			JetColorMap(t, rgb);

			/* 壁粒子は少し暗くして流体と見分けやすくする */
			float dim = (p.type == PT_WALL) ? 0.55f : 1.0f;
			inst.r = rgb[0] * dim;
			inst.g = rgb[1] * dim;
			inst.b = rgb[2] * dim;
		}

		g_instances.push_back(inst);
	}

	return (int)g_instances.size();
}

/*=====================================================================
  【関数名】CreateRenderTargets
  【機能】  スワップチェーンの裏画面 (バックバッファ) から描画先ビューを
            作り、同じ大きさの Z バッファ (深度ステンシル) も作る。
            ウィンドウサイズが変わるたびに作り直す必要がある。
  【引数】  なし (void)
  【戻り値】bool : 成功したら true
  【呼び出し元】InitD3D(), OnResize()
=====================================================================*/
static bool CreateRenderTargets(void)
{
	/* --- 裏画面を取得して「描画先ビュー」を作る --- */
	ID3D11Texture2D* backBuffer = nullptr;
	if (FAILED(g_swapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (void**)&backBuffer))) {
		fprintf(stderr, "ERROR: GetBuffer failed.\n");
		return false;
	}

	D3D11_TEXTURE2D_DESC backBufferDesc;
	backBuffer->GetDesc(&backBufferDesc);
	g_backBufferWidth  = (int)backBufferDesc.Width;
	g_backBufferHeight = (int)backBufferDesc.Height;

	HRESULT hr = g_device->CreateRenderTargetView(backBuffer, nullptr, &g_renderTargetView);
	backBuffer->Release();
	if (FAILED(hr)) {
		fprintf(stderr, "ERROR: CreateRenderTargetView failed.\n");
		return false;
	}

	/* --- Z バッファ (深度バッファ) を作る ---
	   3 次元計算の結果を回して見るとき、奥の粒子が手前に
	   描かれてしまわないようにするために必要 */
	D3D11_TEXTURE2D_DESC depthDesc = {};
	depthDesc.Width      = backBufferDesc.Width;
	depthDesc.Height     = backBufferDesc.Height;
	depthDesc.MipLevels  = 1;
	depthDesc.ArraySize  = 1;
	depthDesc.Format     = DXGI_FORMAT_D32_FLOAT;
	depthDesc.SampleDesc.Count = 1;
	depthDesc.Usage      = D3D11_USAGE_DEFAULT;
	depthDesc.BindFlags  = D3D11_BIND_DEPTH_STENCIL;

	if (FAILED(g_device->CreateTexture2D(&depthDesc, nullptr, &g_depthStencil))) {
		fprintf(stderr, "ERROR: CreateTexture2D (depth) failed.\n");
		return false;
	}
	if (FAILED(g_device->CreateDepthStencilView(g_depthStencil, nullptr, &g_depthStencilView))) {
		fprintf(stderr, "ERROR: CreateDepthStencilView failed.\n");
		return false;
	}
	return true;
}

/*=====================================================================
  【関数名】ReleaseRenderTargets
  【機能】  描画先ビューと Z バッファを解放する。
            ウィンドウサイズ変更時、ResizeBuffers() を呼ぶ前に
            必ずこれらを解放しておかなければならない。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】OnResize(), Cleanup()
=====================================================================*/
static void ReleaseRenderTargets(void)
{
	SAFE_RELEASE(g_depthStencilView);
	SAFE_RELEASE(g_depthStencil);
	SAFE_RELEASE(g_renderTargetView);
}

/*=====================================================================
  【関数名】CompileShaders
  【機能】  埋め込んである HLSL を実行時コンパイルし、
            頂点シェーダ・ピクセルシェーダと入力レイアウトを作る。

            入力レイアウトは「頂点バッファの中身がどう並んでいるか」を
            GPU に教えるもの。ここでは 2 本のバッファを使う。
              スロット 0 : 正方形の隅 (頂点ごとに変わる)
              スロット 1 : 粒子の中心と色 (粒子 1 個ごとに変わる)
            後者に D3D11_INPUT_PER_INSTANCE_DATA を指定することで
            インスタンシング描画ができるようになる。

  【引数】  なし (void)
  【戻り値】bool : 成功したら true
  【呼び出し元】InitD3D()
=====================================================================*/
static bool CompileShaders(void)
{
	UINT flags = D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_OPTIMIZATION_LEVEL3;

	ID3DBlob* vsBlob = nullptr;
	ID3DBlob* psBlob = nullptr;
	ID3DBlob* errBlob = nullptr;

	/* --- 頂点シェーダ --- */
	HRESULT hr = D3DCompile(kShaderSource, strlen(kShaderSource), "mps_shader",
	                        nullptr, nullptr, "VSMain", "vs_5_0", flags, 0,
	                        &vsBlob, &errBlob);
	if (FAILED(hr)) {
		fprintf(stderr, "ERROR: vertex shader compile failed.\n%s\n",
		        errBlob ? (const char*)errBlob->GetBufferPointer() : "(no message)");
		SAFE_RELEASE(errBlob);
		return false;
	}
	SAFE_RELEASE(errBlob);

	/* --- ピクセルシェーダ --- */
	hr = D3DCompile(kShaderSource, strlen(kShaderSource), "mps_shader",
	                nullptr, nullptr, "PSMain", "ps_5_0", flags, 0,
	                &psBlob, &errBlob);
	if (FAILED(hr)) {
		fprintf(stderr, "ERROR: pixel shader compile failed.\n%s\n",
		        errBlob ? (const char*)errBlob->GetBufferPointer() : "(no message)");
		SAFE_RELEASE(errBlob);
		SAFE_RELEASE(vsBlob);
		return false;
	}
	SAFE_RELEASE(errBlob);

	if (FAILED(g_device->CreateVertexShader(vsBlob->GetBufferPointer(),
	                                        vsBlob->GetBufferSize(), nullptr, &g_vertexShader)) ||
	    FAILED(g_device->CreatePixelShader(psBlob->GetBufferPointer(),
	                                       psBlob->GetBufferSize(), nullptr, &g_pixelShader))) {
		fprintf(stderr, "ERROR: CreateVertex/PixelShader failed.\n");
		SAFE_RELEASE(vsBlob); SAFE_RELEASE(psBlob);
		return false;
	}

	/* --- 入力レイアウト --- */
	D3D11_INPUT_ELEMENT_DESC layout[] = {
		/* セマンティクス, 番号, 形式, スロット, オフセット, 種別, ステップ率 */
		{ "POSITION", 0, DXGI_FORMAT_R32G32_FLOAT,       0,  0, D3D11_INPUT_PER_VERTEX_DATA,   0 },
		{ "INSTPOS",  0, DXGI_FORMAT_R32G32B32_FLOAT,    1,  0, D3D11_INPUT_PER_INSTANCE_DATA, 1 },
		{ "INSTCOL",  0, DXGI_FORMAT_R32G32B32A32_FLOAT, 1, 12, D3D11_INPUT_PER_INSTANCE_DATA, 1 },
	};
	hr = g_device->CreateInputLayout(layout, 3,
	                                 vsBlob->GetBufferPointer(), vsBlob->GetBufferSize(),
	                                 &g_inputLayout);
	SAFE_RELEASE(vsBlob);
	SAFE_RELEASE(psBlob);
	if (FAILED(hr)) {
		fprintf(stderr, "ERROR: CreateInputLayout failed.\n");
		return false;
	}
	return true;
}

/*=====================================================================
  【関数名】CreateBuffers
  【機能】  描画に使う 3 本のバッファと、描画設定 (ラスタライザ状態・
            深度ステンシル状態) を作る。

              ・頂点バッファ       : 正方形の 4 隅。内容は不変なので
                                     D3D11_USAGE_IMMUTABLE で作る
              ・インスタンスバッファ: 毎フレーム書き換えるので
                                     D3D11_USAGE_DYNAMIC で作る
              ・定数バッファ       : 行列などを毎フレーム送る

  【引数】  なし (void)
  【戻り値】bool : 成功したら true
  【呼び出し元】InitD3D()
=====================================================================*/
static bool CreateBuffers(void)
{
	/* --- 正方形の 4 隅 (トライアングルストリップの順序) --- */
	const float quadCorners[4][2] = {
		{ -1.0f, -1.0f },
		{ -1.0f,  1.0f },
		{  1.0f, -1.0f },
		{  1.0f,  1.0f },
	};
	D3D11_BUFFER_DESC vbDesc = {};
	vbDesc.ByteWidth = sizeof(quadCorners);
	vbDesc.Usage     = D3D11_USAGE_IMMUTABLE;
	vbDesc.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	D3D11_SUBRESOURCE_DATA vbData = {};
	vbData.pSysMem = quadCorners;
	if (FAILED(g_device->CreateBuffer(&vbDesc, &vbData, &g_quadVertexBuffer))) {
		fprintf(stderr, "ERROR: CreateBuffer (quad) failed.\n");
		return false;
	}

	/* --- インスタンスバッファ (最大粒子数分を確保) --- */
	D3D11_BUFFER_DESC ibDesc = {};
	ibDesc.ByteWidth      = sizeof(InstanceData) * (g_maxParticles > 0 ? g_maxParticles : 1);
	ibDesc.Usage          = D3D11_USAGE_DYNAMIC;
	ibDesc.BindFlags      = D3D11_BIND_VERTEX_BUFFER;
	ibDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	if (FAILED(g_device->CreateBuffer(&ibDesc, nullptr, &g_instanceBuffer))) {
		fprintf(stderr, "ERROR: CreateBuffer (instance) failed.\n");
		return false;
	}

	/* --- 定数バッファ --- */
	D3D11_BUFFER_DESC cbDesc = {};
	cbDesc.ByteWidth      = sizeof(ConstantBufferData);
	cbDesc.Usage          = D3D11_USAGE_DYNAMIC;
	cbDesc.BindFlags      = D3D11_BIND_CONSTANT_BUFFER;
	cbDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	if (FAILED(g_device->CreateBuffer(&cbDesc, nullptr, &g_constantBuffer))) {
		fprintf(stderr, "ERROR: CreateBuffer (constant) failed.\n");
		return false;
	}

	/* --- ラスタライザ状態 ---
	   ビルボードは常にカメラを向くので裏面カリングは不要。
	   カリングを切っておくと頂点の並び順を気にしなくてよい */
	D3D11_RASTERIZER_DESC rsDesc = {};
	rsDesc.FillMode = D3D11_FILL_SOLID;
	rsDesc.CullMode = D3D11_CULL_NONE;
	rsDesc.DepthClipEnable = TRUE;
	if (FAILED(g_device->CreateRasterizerState(&rsDesc, &g_rasterizerState))) {
		fprintf(stderr, "ERROR: CreateRasterizerState failed.\n");
		return false;
	}

	/* --- 深度ステンシル状態 ---
	   粒子は不透明なので、深度テストも深度書き込みも有効にする */
	D3D11_DEPTH_STENCIL_DESC dsDesc = {};
	dsDesc.DepthEnable    = TRUE;
	dsDesc.DepthWriteMask = D3D11_DEPTH_WRITE_MASK_ALL;
	dsDesc.DepthFunc      = D3D11_COMPARISON_LESS;
	if (FAILED(g_device->CreateDepthStencilState(&dsDesc, &g_depthState))) {
		fprintf(stderr, "ERROR: CreateDepthStencilState failed.\n");
		return false;
	}
	return true;
}

/*=====================================================================
  【関数名】InitD3D
  【機能】  DirectX 11 を初期化する。
            デバイス (リソースを作る係) とデバイスコンテキスト
            (描画命令を出す係)、スワップチェーン (画面への出力) を
            まとめて作り、続いて描画先・シェーダ・バッファを用意する。
  【引数】  なし (void)
  【戻り値】bool : 成功したら true
  【呼び出し元】main()
=====================================================================*/
static bool InitD3D(void)
{
	/* --- スワップチェーンの設定 --- */
	DXGI_SWAP_CHAIN_DESC scDesc = {};
	scDesc.BufferCount       = 2;                                  /* 裏画面の枚数 */
	scDesc.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;         /* RGBA 各 8bit */
	scDesc.BufferUsage       = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	scDesc.OutputWindow      = g_hWnd;
	scDesc.SampleDesc.Count  = 1;                                  /* MSAA なし    */
	scDesc.Windowed          = TRUE;
	scDesc.SwapEffect        = DXGI_SWAP_EFFECT_DISCARD;
	/* 幅・高さを 0 にしておくとウィンドウのクライアント領域に合わせてくれる */

	/* ハードウェア (GPU) を優先し、だめならソフトウェア実装 (WARP) で動かす。
	   WARP は遅いが、GPU が無い環境でも確実に描画できる */
	const D3D_DRIVER_TYPE driverTypes[] = {
		D3D_DRIVER_TYPE_HARDWARE,
		D3D_DRIVER_TYPE_WARP,
	};
	const D3D_FEATURE_LEVEL featureLevels[] = {
		D3D_FEATURE_LEVEL_11_1,
		D3D_FEATURE_LEVEL_11_0,
	};

	HRESULT hr = E_FAIL;
	D3D_FEATURE_LEVEL obtainedLevel = D3D_FEATURE_LEVEL_11_0;

	for (int i = 0; i < 2; i++) {
		hr = D3D11CreateDeviceAndSwapChain(
			nullptr,                    /* 既定のアダプタを使う            */
			driverTypes[i],
			nullptr,                    /* ソフトウェアラスタライザ無し    */
			0,                          /* フラグ (デバッグ層は使わない)   */
			featureLevels, 2,
			D3D11_SDK_VERSION,
			&scDesc,
			&g_swapChain, &g_device, &obtainedLevel, &g_context);
		if (SUCCEEDED(hr)) {
			printf("  Direct3D 11 device created (%s, feature level %s)\n",
			       (i == 0) ? "HARDWARE" : "WARP (software)",
			       (obtainedLevel == D3D_FEATURE_LEVEL_11_1) ? "11.1" : "11.0");
			break;
		}
	}
	if (FAILED(hr)) {
		fprintf(stderr, "ERROR: D3D11CreateDeviceAndSwapChain failed (hr=0x%08lX).\n", (unsigned long)hr);
		return false;
	}

	/* Alt+Enter による勝手な全画面切り替えを無効にする */
	{
		IDXGIDevice*  dxgiDevice  = nullptr;
		IDXGIAdapter* dxgiAdapter = nullptr;
		IDXGIFactory* dxgiFactory = nullptr;
		if (SUCCEEDED(g_device->QueryInterface(__uuidof(IDXGIDevice), (void**)&dxgiDevice))) {
			if (SUCCEEDED(dxgiDevice->GetAdapter(&dxgiAdapter))) {
				if (SUCCEEDED(dxgiAdapter->GetParent(__uuidof(IDXGIFactory), (void**)&dxgiFactory))) {
					dxgiFactory->MakeWindowAssociation(g_hWnd, DXGI_MWA_NO_ALT_ENTER);
					dxgiFactory->Release();
				}
				dxgiAdapter->Release();
			}
			dxgiDevice->Release();
		}
	}

	if (!CreateRenderTargets()) return false;
	if (!CompileShaders())      return false;
	if (!CreateBuffers())       return false;
	return true;
}

/*=====================================================================
  【関数名】OnResize
  【機能】  ウィンドウサイズが変わったときに、裏画面と Z バッファを
            新しいサイズで作り直す。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】WndProc()  (WM_SIZE メッセージ受信時)
=====================================================================*/
static void OnResize(void)
{
	if (g_swapChain == nullptr) return;

	/* 裏画面を参照しているビューを先に全部手放す必要がある */
	g_context->OMSetRenderTargets(0, nullptr, nullptr);
	ReleaseRenderTargets();

	/* 幅・高さに 0 を渡すと現在のウィンドウサイズに合わせてくれる */
	g_swapChain->ResizeBuffers(0, 0, 0, DXGI_FORMAT_UNKNOWN, 0);

	CreateRenderTargets();
}

/*=====================================================================
  【関数名】RenderFrame
  【機能】  1 コマ分の描画を行う。

            (1) カメラの位置と向きからビュー行列・射影行列を作る
            (2) 定数バッファを更新して GPU に送る
            (3) 現在フレームの粒子からインスタンスデータを作って送る
            (4) 画面を塗りつぶし、DrawInstanced() で全粒子を一気に描く

  【引数】  bool present : true なら描画後に画面へ表示する。
                           スクリーンショット用に false で呼ぶと
                           裏画面に描くだけで表示しない。
  【戻り値】なし (void)
  【呼び出し元】MainLoop(), SaveScreenshotMode()
=====================================================================*/
static void RenderFrame(bool present)
{
	if (g_renderTargetView == nullptr) return;

	/*--- (1) カメラ ---------------------------------------------*/

	/* 視線方向 (前方ベクトル)。ヨー角とピッチ角から求める */
	XMVECTOR forward = XMVectorSet(
		cosf(g_cameraPitch) * sinf(g_cameraYaw),
		sinf(g_cameraPitch),
		cosf(g_cameraPitch) * cosf(g_cameraYaw),
		0.0f);

	XMVECTOR target = XMLoadFloat3(&g_cameraTarget);
	XMVECTOR eye    = XMVectorSubtract(target, XMVectorScale(forward, g_cameraDistance));

	/* 左手座標系での 右 = 上 × 前 、 上 = 前 × 右 */
	XMVECTOR worldUp = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);
	XMVECTOR right   = XMVector3Normalize(XMVector3Cross(worldUp, forward));
	XMVECTOR up      = XMVector3Cross(forward, right);

	XMMATRIX view = XMMatrixLookAtLH(eye, target, up);

	/* 正射影 (遠近感を付けない)。
	   2 次元計算の結果を見るときは遠近法が無いほうが形が正しく見える */
	float aspect = (g_backBufferHeight > 0)
	             ? (float)g_backBufferWidth / (float)g_backBufferHeight : 1.0f;
	XMMATRIX proj = XMMatrixOrthographicLH(g_orthoHeight * aspect, g_orthoHeight,
	                                       0.01f, g_cameraDistance * 4.0f + 100.0f);

	/*--- (2) 定数バッファの更新 ----------------------------------*/
	{
		D3D11_MAPPED_SUBRESOURCE mapped;
		if (SUCCEEDED(g_context->Map(g_constantBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
			ConstantBufferData* cb = (ConstantBufferData*)mapped.pData;

			/* HLSL の行列は既定で列優先なので、転置して渡す */
			XMStoreFloat4x4(&cb->viewProj, XMMatrixTranspose(XMMatrixMultiply(view, proj)));

			/* ビルボードを広げる量 = カメラの右/上方向 × 粒子半径 */
			XMFLOAT3 r3, u3;
			XMStoreFloat3(&r3, XMVectorScale(right, g_particleRadius));
			XMStoreFloat3(&u3, XMVectorScale(up,    g_particleRadius));
			cb->right = XMFLOAT4(r3.x, r3.y, r3.z, 0.0f);
			cb->up    = XMFLOAT4(u3.x, u3.y, u3.z, 0.0f);

			g_context->Unmap(g_constantBuffer, 0);
		}
	}

	/*--- (3) インスタンスデータの更新 ----------------------------*/
	int frameIndex = (int)g_currentFrame;
	if (frameIndex < 0) frameIndex = 0;
	if (frameIndex >= (int)g_frames.size()) frameIndex = (int)g_frames.size() - 1;

	int instanceCount = BuildInstances(frameIndex);
	if (instanceCount > 0) {
		D3D11_MAPPED_SUBRESOURCE mapped;
		if (SUCCEEDED(g_context->Map(g_instanceBuffer, 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped))) {
			memcpy(mapped.pData, g_instances.data(), sizeof(InstanceData) * instanceCount);
			g_context->Unmap(g_instanceBuffer, 0);
		}
	}

	/*--- (4) 描画 ------------------------------------------------*/

	/* 描画先とビューポートを設定する */
	g_context->OMSetRenderTargets(1, &g_renderTargetView, g_depthStencilView);

	D3D11_VIEWPORT viewport = {};
	viewport.Width    = (float)g_backBufferWidth;
	viewport.Height   = (float)g_backBufferHeight;
	viewport.MaxDepth = 1.0f;
	g_context->RSSetViewports(1, &viewport);

	/* 画面を暗い紺色で塗りつぶし、Z バッファを 1.0 で初期化する */
	const float clearColor[4] = { 0.07f, 0.08f, 0.12f, 1.0f };
	g_context->ClearRenderTargetView(g_renderTargetView, clearColor);
	if (g_depthStencilView) {
		g_context->ClearDepthStencilView(g_depthStencilView, D3D11_CLEAR_DEPTH, 1.0f, 0);
	}

	if (instanceCount > 0) {
		/* 2 本の頂点バッファを同時に設定する。
		   スロット 0 = 正方形の隅、スロット 1 = 粒子ごとの情報 */
		ID3D11Buffer* buffers[2] = { g_quadVertexBuffer, g_instanceBuffer };
		UINT strides[2] = { sizeof(float) * 2, sizeof(InstanceData) };
		UINT offsets[2] = { 0, 0 };
		g_context->IASetVertexBuffers(0, 2, buffers, strides, offsets);
		g_context->IASetInputLayout(g_inputLayout);
		g_context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);

		g_context->VSSetShader(g_vertexShader, nullptr, 0);
		g_context->VSSetConstantBuffers(0, 1, &g_constantBuffer);
		g_context->PSSetShader(g_pixelShader, nullptr, 0);

		g_context->RSSetState(g_rasterizerState);
		g_context->OMSetDepthStencilState(g_depthState, 0);

		/* 4 頂点の正方形を instanceCount 個ぶん描く = 全粒子を 1 命令で描画 */
		g_context->DrawInstanced(4, instanceCount, 0, 0);
	}

	if (present) {
		g_swapChain->Present(1, 0);   /* 1 = 垂直同期あり (ティアリング防止) */
	}
}

/*=====================================================================
  【関数名】SaveBackBufferAsBMP
  【機能】  裏画面 (バックバッファ) の内容を 24bit BMP として保存する。

            GPU 上のテクスチャは CPU から直接読めないので、
            CPU が読める「ステージングテクスチャ」を作って
            CopyResource() で内容をコピーしてから読み出す。

            BMP は行が下から上へ並ぶ形式なので、逆順に書き出す。

  【引数】  const char* path : 保存先のファイルパス
  【戻り値】bool : 保存できたら true
  【呼び出し元】WndProc() (F12 キー), SaveScreenshotMode()
=====================================================================*/
static bool SaveBackBufferAsBMP(const char* path)
{
	ID3D11Texture2D* backBuffer = nullptr;
	if (FAILED(g_swapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), (void**)&backBuffer))) {
		return false;
	}

	D3D11_TEXTURE2D_DESC desc;
	backBuffer->GetDesc(&desc);

	/* CPU が読み出せるテクスチャを用意する */
	D3D11_TEXTURE2D_DESC stagingDesc = desc;
	stagingDesc.Usage          = D3D11_USAGE_STAGING;
	stagingDesc.BindFlags      = 0;
	stagingDesc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;
	stagingDesc.MiscFlags      = 0;

	ID3D11Texture2D* staging = nullptr;
	if (FAILED(g_device->CreateTexture2D(&stagingDesc, nullptr, &staging))) {
		backBuffer->Release();
		return false;
	}
	g_context->CopyResource(staging, backBuffer);
	backBuffer->Release();

	D3D11_MAPPED_SUBRESOURCE mapped;
	if (FAILED(g_context->Map(staging, 0, D3D11_MAP_READ, 0, &mapped))) {
		staging->Release();
		return false;
	}

	const int width  = (int)desc.Width;
	const int height = (int)desc.Height;
	const int rowBytes = width * 3;
	const int padding  = (4 - (rowBytes % 4)) % 4;   /* BMP の行は 4 バイト境界 */
	const int rowSize  = rowBytes + padding;

	FILE* fp = fopen(path, "wb");
	if (fp == nullptr) {
		g_context->Unmap(staging, 0);
		staging->Release();
		return false;
	}

	BITMAPFILEHEADER fileHeader = {};
	BITMAPINFOHEADER infoHeader = {};
	fileHeader.bfType    = 0x4D42;   /* 'BM' */
	fileHeader.bfOffBits = sizeof(BITMAPFILEHEADER) + sizeof(BITMAPINFOHEADER);
	fileHeader.bfSize    = fileHeader.bfOffBits + rowSize * height;
	infoHeader.biSize     = sizeof(BITMAPINFOHEADER);
	infoHeader.biWidth    = width;
	infoHeader.biHeight   = height;     /* 正の値 = 下から上へ並ぶ */
	infoHeader.biPlanes   = 1;
	infoHeader.biBitCount = 24;
	infoHeader.biCompression = BI_RGB;

	fwrite(&fileHeader, sizeof(fileHeader), 1, fp);
	fwrite(&infoHeader, sizeof(infoHeader), 1, fp);

	std::vector<unsigned char> row(rowSize, 0);
	const unsigned char* src = (const unsigned char*)mapped.pData;

	for (int y = height - 1; y >= 0; y--) {   /* BMP は最下行から書く */
		const unsigned char* line = src + (size_t)y * mapped.RowPitch;
		for (int x = 0; x < width; x++) {
			/* 裏画面は RGBA の順、BMP は BGR の順 */
			row[x * 3 + 0] = line[x * 4 + 2];  /* B */
			row[x * 3 + 1] = line[x * 4 + 1];  /* G */
			row[x * 3 + 2] = line[x * 4 + 0];  /* R */
		}
		fwrite(row.data(), 1, rowSize, fp);
	}

	fclose(fp);
	g_context->Unmap(staging, 0);
	staging->Release();
	return true;
}

/*=====================================================================
  【関数名】UpdateWindowTitle
  【機能】  現在の再生状態をウィンドウのタイトルバーに表示する。
            画面に文字を描く代わりの簡易的な情報表示。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】MainLoop()  (約 0.1 秒ごと)
=====================================================================*/
static void UpdateWindowTitle(void)
{
	int frameIndex = (int)g_currentFrame;
	if (frameIndex < 0) frameIndex = 0;
	if (frameIndex >= (int)g_frames.size()) frameIndex = (int)g_frames.size() - 1;

	char title[512];
	sprintf(title,
	        "MPS Viewer (DirectX 11)  |  Frame %d/%d  t=%.3f s  |  Color: %s (x%.2f)  |  Speed x%.2f  |  %s%s",
	        frameIndex, (int)g_frames.size() - 1,
	        g_frames[frameIndex].time,
	        kColorModeName[g_colorMode], g_colorScale,
	        g_playbackSpeed,
	        g_isPlaying ? "Playing" : "Paused",
	        g_isLooping ? " [Loop]" : "");
	SetWindowTextA(g_hWnd, title);
}

/*=====================================================================
  【関数名】WndProc
  【機能】  ウィンドウプロシージャ。Windows から送られてくる
            メッセージ (キー入力・マウス操作・サイズ変更など) を
            処理する。Win32 プログラムの中心となる関数。
  【引数】  HWND   hWnd    : メッセージの宛先ウィンドウ
            UINT   message : メッセージの種類 (WM_KEYDOWN など)
            WPARAM wParam  : メッセージごとの付加情報 (キーコードなど)
            LPARAM lParam  : メッセージごとの付加情報 (座標など)
  【戻り値】LRESULT : 処理結果。自分で処理しないものは
                      DefWindowProc() に任せてその戻り値を返す。
  【呼び出し元】Windows (DispatchMessage() 経由で呼ばれる)
=====================================================================*/
static LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	switch (message) {

	case WM_DESTROY:
		PostQuitMessage(0);      /* メッセージループを終わらせる */
		return 0;

	case WM_SIZE:
		if (g_device && wParam != SIZE_MINIMIZED) OnResize();
		return 0;

	/*--- キーボード ---------------------------------------------*/
	case WM_KEYDOWN:
		switch (wParam) {
		case VK_ESCAPE:
			PostMessage(hWnd, WM_CLOSE, 0, 0);
			break;

		case VK_SPACE:
			g_isPlaying = !g_isPlaying;
			break;

		case VK_RIGHT:   /* 1 コマ進む */
			g_isPlaying = false;
			g_currentFrame = (float)((int)g_currentFrame + 1);
			if (g_currentFrame > (float)(g_frames.size() - 1)) {
				g_currentFrame = (float)(g_frames.size() - 1);
			}
			break;

		case VK_LEFT:    /* 1 コマ戻る */
			g_isPlaying = false;
			g_currentFrame = (float)((int)g_currentFrame - 1);
			if (g_currentFrame < 0.0f) g_currentFrame = 0.0f;
			break;

		case VK_HOME:
			g_currentFrame = 0.0f;
			break;

		case VK_F12:
			{
				char name[256];
				sprintf(name, "shot_%03d.bmp", g_screenshotCount);
				if (SaveBackBufferAsBMP(name)) {
					printf("  screenshot saved: %s\n", name);
					g_screenshotCount++;
				} else {
					printf("  screenshot FAILED\n");
				}
			}
			break;

		case 'L':
			g_isLooping = !g_isLooping;
			break;

		case 'W':
			g_showWall = !g_showWall;
			break;

		case 'D':
			g_showDummyWall = !g_showDummyWall;
			break;

		case 'F':
			FitCameraToBounds();
			break;

		case '1': g_colorMode = CM_PRESSURE; break;
		case '2': g_colorMode = CM_VELOCITY; break;
		case '3': g_colorMode = CM_TYPE;     break;
		case '4': g_colorMode = CM_DENSITY;  break;

		case '0':
			g_playbackSpeed = 1.0f;
			break;

		case VK_OEM_PLUS:  case VK_ADD:
			g_playbackSpeed = Clampf(g_playbackSpeed * 1.25f, 0.05f, 20.0f);
			break;

		case VK_OEM_MINUS: case VK_SUBTRACT:
			g_playbackSpeed = Clampf(g_playbackSpeed / 1.25f, 0.05f, 20.0f);
			break;

		case VK_OEM_4:   /* [ キー : カラースケール上限を下げる */
			g_colorScale = Clampf(g_colorScale * 0.8f, 0.01f, 100.0f);
			break;

		case VK_OEM_6:   /* ] キー : カラースケール上限を上げる */
			g_colorScale = Clampf(g_colorScale * 1.25f, 0.01f, 100.0f);
			break;

		case VK_OEM_COMMA:   /* , キー : 粒子を小さく */
			g_particleRadius = Clampf(g_particleRadius * 0.9f,
			                          g_particleSpacing * 0.05f, g_particleSpacing * 3.0f);
			break;

		case VK_OEM_PERIOD:  /* . キー : 粒子を大きく */
			g_particleRadius = Clampf(g_particleRadius * 1.1f,
			                          g_particleSpacing * 0.05f, g_particleSpacing * 3.0f);
			break;
		}
		return 0;

	/*--- マウス -------------------------------------------------*/
	case WM_LBUTTONDOWN:
		g_isDraggingLeft = true;
		g_lastMousePos.x = GET_X_LPARAM_COMPAT(lParam);
		g_lastMousePos.y = GET_Y_LPARAM_COMPAT(lParam);
		SetCapture(hWnd);
		return 0;

	case WM_RBUTTONDOWN:
		g_isDraggingRight = true;
		g_lastMousePos.x = GET_X_LPARAM_COMPAT(lParam);
		g_lastMousePos.y = GET_Y_LPARAM_COMPAT(lParam);
		SetCapture(hWnd);
		return 0;

	case WM_LBUTTONUP:
		g_isDraggingLeft = false;
		if (!g_isDraggingRight) ReleaseCapture();
		return 0;

	case WM_RBUTTONUP:
		g_isDraggingRight = false;
		if (!g_isDraggingLeft) ReleaseCapture();
		return 0;

	case WM_MOUSEMOVE:
		{
			int mouseX = GET_X_LPARAM_COMPAT(lParam);
			int mouseY = GET_Y_LPARAM_COMPAT(lParam);
			int dx = mouseX - g_lastMousePos.x;
			int dy = mouseY - g_lastMousePos.y;
			g_lastMousePos.x = mouseX;
			g_lastMousePos.y = mouseY;

			if (g_isDraggingLeft) {
				/* 回転 : 横移動でヨー角、縦移動でピッチ角を変える */
				g_cameraYaw   -= dx * 0.008f;
				g_cameraPitch = Clampf(g_cameraPitch + dy * 0.008f, -1.50f, 1.50f);
			}
			if (g_isDraggingRight) {
				/* 平行移動 : 画面上の移動量を world 単位に換算する */
				float unitsPerPixel = g_orthoHeight / (float)(g_backBufferHeight > 0 ? g_backBufferHeight : 1);

				XMVECTOR forward = XMVectorSet(
					cosf(g_cameraPitch) * sinf(g_cameraYaw),
					sinf(g_cameraPitch),
					cosf(g_cameraPitch) * cosf(g_cameraYaw), 0.0f);
				XMVECTOR worldUp = XMVectorSet(0.0f, 1.0f, 0.0f, 0.0f);
				XMVECTOR right   = XMVector3Normalize(XMVector3Cross(worldUp, forward));
				XMVECTOR up      = XMVector3Cross(forward, right);

				XMVECTOR target = XMLoadFloat3(&g_cameraTarget);
				target = XMVectorSubtract(target, XMVectorScale(right, dx * unitsPerPixel));
				target = XMVectorAdd(target,      XMVectorScale(up,    dy * unitsPerPixel));
				XMStoreFloat3(&g_cameraTarget, target);
			}
		}
		return 0;

	case WM_MOUSEWHEEL:
		{
			int delta = GET_WHEEL_DELTA_WPARAM(wParam);
			float factor = (delta > 0) ? (1.0f / 1.15f) : 1.15f;  /* 手前で拡大 */
			g_orthoHeight = Clampf(g_orthoHeight * factor, 1e-4f, 1e6f);
		}
		return 0;
	}

	return DefWindowProc(hWnd, message, wParam, lParam);
}

/*=====================================================================
  【関数名】CreateAppWindow
  【機能】  描画用の Win32 ウィンドウを作る。
            ウィンドウクラスを登録し、指定したクライアント領域の
            大きさになるようにウィンドウを作成する。
  【引数】  int width  : クライアント領域の幅   [pixel]
            int height : クライアント領域の高さ [pixel]
            bool visible : 表示するかどうか
                           (スクリーンショット専用モードでは false)
  【戻り値】bool : 成功したら true
  【副作用】g_hWnd にウィンドウハンドルを設定する
  【呼び出し元】main()
=====================================================================*/
static bool CreateAppWindow(int width, int height, bool visible)
{
	HINSTANCE hInstance = GetModuleHandle(nullptr);

	WNDCLASSEXA wc = {};
	wc.cbSize        = sizeof(wc);
	wc.style         = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc   = WndProc;
	wc.hInstance     = hInstance;
	wc.hCursor       = LoadCursor(nullptr, IDC_ARROW);
	wc.hbrBackground = nullptr;      /* 背景は DirectX が塗るので不要 */
	wc.lpszClassName = "MpsViewerWindowClass";
	if (!RegisterClassExA(&wc)) {
		fprintf(stderr, "ERROR: RegisterClassEx failed.\n");
		return false;
	}

	/* 指定したクライアント領域の大きさになるよう、枠の分を足す */
	RECT rect = { 0, 0, width, height };
	AdjustWindowRect(&rect, WS_OVERLAPPEDWINDOW, FALSE);

	g_hWnd = CreateWindowExA(
		0, wc.lpszClassName, "MPS Viewer (DirectX 11)",
		WS_OVERLAPPEDWINDOW,
		CW_USEDEFAULT, CW_USEDEFAULT,
		rect.right - rect.left, rect.bottom - rect.top,
		nullptr, nullptr, hInstance, nullptr);

	if (g_hWnd == nullptr) {
		fprintf(stderr, "ERROR: CreateWindowEx failed.\n");
		return false;
	}

	/* ShowWindow() を 2 回呼んでいるのは Windows の仕様への対処である。

	   プロセスの「最初の」ShowWindow() 呼び出しは、引数で渡した値ではなく
	   プロセス起動時に親から渡された STARTUPINFO の wShowWindow が
	   使われることがある (STARTF_USESHOWWINDOW が立っている場合)。
	   たとえば PowerShell の
	       Start-Process mps_viewer.exe -WindowStyle Hidden
	   のようにバッチやスクリプトから非表示指定で起動されると、
	   こちらが SW_SHOW を渡しても最初の 1 回は無視されて
	   ウィンドウが表示されない。
	   2 回目以降の呼び出しは引数どおりに動くので、もう一度呼んでおく。 */
	const int showCommand = visible ? SW_SHOW : SW_HIDE;
	ShowWindow(g_hWnd, showCommand);
	ShowWindow(g_hWnd, showCommand);

	UpdateWindow(g_hWnd);
	return true;
}

/*=====================================================================
  【関数名】UpdateAnimation
  【機能】  前回の描画からの経過時間をもとにフレーム番号を進める。

            .prof は計算時間 g_frameDeltaTime 秒おきに出力されている
            ので、実時間 dt 秒が経過したら
              dt × 再生速度 ÷ g_frameDeltaTime
            コマだけ進めると、シミュレーションが実時間どおりの速さで
            再生される。

  【引数】  float deltaSeconds : 前回の描画からの実経過時間 [s]
  【戻り値】なし (void)
  【副作用】g_currentFrame を進める
  【呼び出し元】MainLoop()
=====================================================================*/
static void UpdateAnimation(float deltaSeconds)
{
	if (!g_isPlaying)        return;
	if (g_frames.size() < 2) return;

	g_currentFrame += (deltaSeconds * g_playbackSpeed) / g_frameDeltaTime;

	const float lastFrame = (float)(g_frames.size() - 1);
	if (g_currentFrame > lastFrame) {
		if (g_isLooping) {
			g_currentFrame = 0.0f;             /* 先頭に戻って繰り返す */
		} else {
			g_currentFrame = lastFrame;
			g_isPlaying = false;               /* 末尾で停止           */
		}
	}
}

/*=====================================================================
  【関数名】MainLoop
  【機能】  メッセージループ。ウィンドウが閉じられるまで
            「メッセージ処理 → アニメーション更新 → 描画」を繰り返す。

            GetMessage() ではなく PeekMessage() を使うのは、
            メッセージが無いときでも描画を続けるため
            (ゲームやビューアで使われる定番の書き方)。

  【引数】  なし (void)
  【戻り値】int : 終了コード
  【呼び出し元】main()
=====================================================================*/
static int MainLoop(void)
{
	/* 高分解能タイマーの準備 (経過時間を正確に測るため) */
	LARGE_INTEGER frequency, previousCount, currentCount;
	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&previousCount);

	double titleTimer = 0.0;   /* タイトル更新の間隔を測る */

	MSG msg = {};
	while (msg.message != WM_QUIT) {

		/* たまっているメッセージを全部処理する */
		if (PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)) {
			TranslateMessage(&msg);
			DispatchMessage(&msg);
			continue;
		}

		/* 前回からの経過時間を求める */
		QueryPerformanceCounter(&currentCount);
		double deltaSeconds = (double)(currentCount.QuadPart - previousCount.QuadPart)
		                    / (double)frequency.QuadPart;
		previousCount = currentCount;

		/* ウィンドウを動かしている間など、極端に大きい値は切り捨てる */
		if (deltaSeconds > 0.25) deltaSeconds = 0.25;

		UpdateAnimation((float)deltaSeconds);
		RenderFrame(true);

		/* タイトルバーは 0.1 秒に 1 回だけ更新する (毎回だと重い) */
		titleTimer += deltaSeconds;
		if (titleTimer > 0.1) {
			UpdateWindowTitle();
			titleTimer = 0.0;
		}
	}
	return (int)msg.wParam;
}

/*=====================================================================
  【関数名】Cleanup
  【機能】  DirectX のオブジェクトをすべて解放する。
            作ったのと逆の順序で解放するのが基本。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】main()
=====================================================================*/
static void Cleanup(void)
{
	if (g_context) g_context->ClearState();

	SAFE_RELEASE(g_depthState);
	SAFE_RELEASE(g_rasterizerState);
	SAFE_RELEASE(g_constantBuffer);
	SAFE_RELEASE(g_instanceBuffer);
	SAFE_RELEASE(g_quadVertexBuffer);
	SAFE_RELEASE(g_inputLayout);
	SAFE_RELEASE(g_pixelShader);
	SAFE_RELEASE(g_vertexShader);
	ReleaseRenderTargets();
	SAFE_RELEASE(g_swapChain);
	SAFE_RELEASE(g_context);
	SAFE_RELEASE(g_device);
}

/*=====================================================================
  【関数名】PrintUsage
  【機能】  コンソールに操作説明を表示する。
  【引数】  なし (void)
  【戻り値】なし (void)
  【呼び出し元】main()
=====================================================================*/
static void PrintUsage(void)
{
	printf("\n");
	printf("  ---------------- KEY / MOUSE ----------------\n");
	printf("   Space      : play / pause\n");
	printf("   Right/Left : step one frame\n");
	printf("   Home       : jump to first frame\n");
	printf("   L          : loop on/off\n");
	printf("   + / -      : playback speed up / down     0 : reset speed\n");
	printf("   1/2/3/4    : color by Pressure / Velocity / Type / NumberDensity\n");
	printf("   [ / ]      : color scale max  down / up\n");
	printf("   W / D      : show wall / dummy-wall particles\n");
	printf("   , / .      : particle radius  smaller / larger\n");
	printf("   F          : fit view to the whole domain\n");
	printf("   F12        : save screenshot (shot_NNN.bmp)\n");
	printf("   Esc        : quit\n");
	printf("   Mouse L-drag : rotate   R-drag : pan   Wheel : zoom\n");
	printf("  ---------------------------------------------\n\n");
}

/*=====================================================================
  【関数名】main
  【機能】  プログラムの入口。
            引数を解釈し、データを読み込み、ウィンドウと DirectX を
            用意してアニメーション再生を始める。

            なお本プログラムはコンソールアプリケーションとして作り、
            ウィンドウは CreateWindowEx() で自分で作っている。
            こうすると printf() による進捗表示やエラー表示が
            そのまま使えて、学習用途では扱いやすい。

  【引数】  int argc    : コマンドライン引数の個数
            char** argv : コマンドライン引数の文字列配列
  【戻り値】int : 正常終了なら 0、エラーなら 1
=====================================================================*/
int main(int argc, char** argv)
{
	printf("\n*** MPS VIEWER (DirectX 11) ***\n");

	/*--- コマンドライン引数の解釈 --------------------------------*/
	const char* directory   = ".";     /* .prof があるフォルダ     */
	int         shotFrame   = -1;      /* -shot で指定したフレーム */
	const char* shotPath    = nullptr; /* -shot の出力ファイル     */
	float       initialYawDeg   = 0.0f;/* -yaw   で指定した水平回転角 [度] */
	float       initialPitchDeg = 0.0f;/* -pitch で指定した仰角       [度] */
	int         wallOption      = -1;  /* 壁の表示指定
                                          -1:指定なし（自動）
                                           0:-nowall で非表示
                                           1:-wall   で表示          */

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-r") == 0 && i + 1 < argc) {
			g_particleRadius = (float)atof(argv[++i]);
		} else if (strcmp(argv[i], "-c") == 0 && i + 1 < argc) {
			int mode = atoi(argv[++i]);
			if (mode >= 0 && mode < CM_COUNT) g_colorMode = mode;
		} else if (strcmp(argv[i], "-wall") == 0) {
			wallOption = 1;
		} else if (strcmp(argv[i], "-nowall") == 0) {
			wallOption = 0;
		} else if (strcmp(argv[i], "-yaw") == 0 && i + 1 < argc) {
			initialYawDeg = (float)atof(argv[++i]);
		} else if (strcmp(argv[i], "-pitch") == 0 && i + 1 < argc) {
			initialPitchDeg = (float)atof(argv[++i]);
		} else if (strcmp(argv[i], "-shot") == 0 && i + 2 < argc) {
			shotFrame = atoi(argv[++i]);
			shotPath  = argv[++i];
		} else if (argv[i][0] != '-') {
			directory = argv[i];
		}
	}

	/*--- データ読み込み ------------------------------------------*/
	printf("  data folder: %s\n", directory);
	if (!LoadAllFrames(directory)) {
		fprintf(stderr,
		        "ERROR: no output_*.prof found in \"%s\".\n"
		        "       Run mps.exe first, or pass the folder as an argument:\n"
		        "         mps_viewer.exe <folder>\n", directory);
		return 1;
	}

	AnalyzeFrames();

	/* --- 壁の表示・非表示を決める ---
	   3 次元計算では壁が閉じた箱になっているため、既定のまま壁を
	   表示すると外側の壁しか見えず、中の流体がまったく確認できない。
	   そこで 3 次元と判定されたときは壁を既定で非表示にする。
	   -wall / -nowall が指定されていればそちらを優先する。         */
	if (wallOption == 0)      g_showWall = false;
	else if (wallOption == 1) g_showWall = true;
	else                      g_showWall = !g_isThreeDimensional;

	printf("  dimension     : %s\n", g_isThreeDimensional ? "3D" : "2D");
	printf("  wall display  : %s%s\n",
	       g_showWall ? "ON" : "OFF",
	       (wallOption < 0 && g_isThreeDimensional)
	           ? "  (hidden automatically for 3D; press W to show)" : "");
	printf("  frames        : %d  (dt = %.4f s per frame)\n",
	       (int)g_frames.size(), g_frameDeltaTime);
	printf("  particles     : %d (max)\n", g_maxParticles);
	printf("  spacing l0    : %.5f m   -> radius %.5f m\n",
	       g_particleSpacing, g_particleRadius);
	printf("  bounds        : (%.3f, %.3f, %.3f) - (%.3f, %.3f, %.3f)\n",
	       g_boundsMin.x, g_boundsMin.y, g_boundsMin.z,
	       g_boundsMax.x, g_boundsMax.y, g_boundsMax.z);
	printf("  pressure      : %.3f .. %.3f Pa   (color scale, 99th percentile)\n",
	       g_pressureMin, g_pressureMax);
	printf("  velocity      : %.3f .. %.3f m/s  (color scale, 99th percentile)\n",
	       g_speedMin, g_speedMax);

	/*--- ウィンドウと DirectX の初期化 ---------------------------*/
	SetProcessDPIAware();   /* 高 DPI 環境でぼやけないようにする */

	const bool isShotMode = (shotPath != nullptr);

	if (!CreateAppWindow(1280, 720, !isShotMode)) return 1;
	if (!InitD3D()) { Cleanup(); return 1; }

	/* 画面の縦横比が確定してからカメラを合わせる
	   (FitCameraToBounds() は縦横比を使うため、D3D 初期化後に呼ぶ) */
	FitCameraToBounds();

	/* -yaw / -pitch が指定されていれば、その向きから見る。
	   3 次元計算の結果を斜めから撮影したいときに使う         */
	{
		const float kDegreeToRadian = 3.14159265f / 180.0f;
		g_cameraYaw   = initialYawDeg * kDegreeToRadian;
		g_cameraPitch = Clampf(initialPitchDeg * kDegreeToRadian, -1.50f, 1.50f);
	}

	/*--- スクリーンショット専用モード ----------------------------*/
	if (isShotMode) {
		if (shotFrame < 0) shotFrame = 0;
		if (shotFrame > (int)g_frames.size() - 1) shotFrame = (int)g_frames.size() - 1;
		g_currentFrame = (float)shotFrame;

		RenderFrame(false);   /* 裏画面に描くだけで画面には出さない */

		bool ok = SaveBackBufferAsBMP(shotPath);
		printf(ok ? "  screenshot saved: %s\n" : "  screenshot FAILED: %s\n", shotPath);

		Cleanup();
		DestroyWindow(g_hWnd);
		return ok ? 0 : 1;
	}

	/*--- 通常のアニメーション再生 --------------------------------*/
	PrintUsage();
	UpdateWindowTitle();

	int exitCode = MainLoop();

	Cleanup();
	printf("*** END ***\n\n");
	return exitCode;
}
