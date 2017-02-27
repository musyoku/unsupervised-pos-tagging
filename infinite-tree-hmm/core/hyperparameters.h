#ifndef _hyper_
#define _hyper_

// [Tree-Structured Stick Breaking for Hierarchical Data](https://hips.seas.harvard.edu/files/adams-tssb-nips-2010.pdf)
// ではこれらのハイパーパラメータは実行前に一度だけ一様分布からサンプリングしていたのでminとmaxを設定する
#define iTHMM_ALPHA_MIN 0.05
#define iTHMM_ALPHA_MAX 0.2
#define iTHMM_GAMMA_MIN 0.05
#define iTHMM_GAMMA_MAX 0.2
#define iTHMM_LAMBDA_MIN 0.05	// 0以上1以下
#define iTHMM_LAMBDA_MAX 0.8	// 0以上1以下

// [無限木構造隠れMarkovモデルによる階層的品詞の教師なし学習](http://chasen.org/~daiti-m/paper/nl226ithmm.pdf)
// ではこれらのハイパーパラメータは固定値
#define iTHMM_TAU_0 1.0
#define iTHMM_TAU_1 100.0

// 以下はHPYLMのハイパーパラメータの初期値
// 実行中にサンプリングして正しい値を推定する
#define HPYLM_D 0.2			// 0以上1未満. ディスカウント係数
#define HPYLM_THETA 2.0		// 集中度
// 以下は上のハイパーパラメータ2つを推定する時に使うハイパーパラメータ
// 実行中常に固定値
#define HPYLM_A 	1.0
#define HPYLM_B 	1.0
#define HPYLM_ALPHA 1.0
#define HPYLM_BETA  1.0

#define EPS 1e-12			// 本来これは0が望ましい

#endif