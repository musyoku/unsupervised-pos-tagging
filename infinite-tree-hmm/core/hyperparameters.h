#ifndef _hyper_
#define _hyper_

// これらのハイパーパラメータは実行前に一度だけ一様分布からサンプリングしていたのでminとmaxを設定する*1
#define iTHMM_ALPHA_MIN 0.1
#define iTHMM_ALPHA_MAX 0.5
#define iTHMM_GAMMA_MIN 0.1
#define iTHMM_GAMMA_MAX 0.5
#define iTHMM_LAMBDA_MIN 0.01	// 0以上1以下
#define iTHMM_LAMBDA_MAX 0.1	// 0以上1以下
// HTSSBで親の情報をどの程度受け継ぐかを制御するパラメータ
// 論文では上のαと同じ記号が使われているが区別する必要がある
#define iTHMM_STRENGTH_MIN 1.0
#define iTHMM_STRENGTH_MAX 1000.0

// これらのハイパーパラメータは固定*2
#define iTHMM_TAU_0 1.0
#define iTHMM_TAU_1 100.0

// 以下はHPYLMのハイパーパラメータの初期値*3
// 実行中にサンプリングして正しい値を推定する
#define HPYLM_D 0.2				// 0以上1未満. ディスカウント係数
#define HPYLM_THETA 0.1			// 集中度
// ルートノードでは単語の出力分布はフラットになってほしいので以下の固定値を用いる
#define HPYLM_D_ROOT 0.2		// ルートノードでのディスカウント係数
#define HPYLM_THETA_ROOT 0.1	// ルートノードでの集中度
// 以下は上のハイパーパラメータ2つを推定する時に使うハイパーパラメータ
// 実行中常に固定値
#define HPYLM_A 	1.0
#define HPYLM_B 	1.0
#define HPYLM_ALPHA 1.0
#define HPYLM_BETA  1.0

#define EPS 1e-12			// 本来これは0が望ましい

// *1 [Tree-Structured Stick Breaking for Hierarchical Data](https://hips.seas.harvard.edu/files/adams-tssb-nips-2010.pdf)
// *2 [無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習](http://chasen.org/~daiti-m/paper/nl226ithmm.pdf)
// *3 [A Bayesian Interpretation of Interpolated Kneser-Ney](https://www.stats.ox.ac.uk/~teh/research/compling/hpylm.pdf)
#endif