/**
 * 出力の各行に対応する．
 * 作物xと作物yを交配してクセSの作物を作ることを表す．
 */
struct cross {
    int x, y;
    char S[52];
};

/**
 * 答えがYESのときyn = 1, NOのときyn = 0で呼び出す．
 * mは交配の回数を表す．
 * cは交配の情報を持ったcross型の配列．
 */
void output(int yn, int m, struct cross *c) {
    if (!yn) {
        printf("NO\n");
    } else {
        printf("YES\n");
        printf("%d\n", m);
        int i;
        for (i = 0; i < m; i++)
            printf("%d %d %s\n", c[i].x, c[i].y, c[i].S);
    }
    fflush(stdout);
}
