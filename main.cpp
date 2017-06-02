#include "bits/stdc++.h"

using namespace std;

/**
 * 出力の各行に対応する．
 * 作物xと作物yを交配してクセSの作物を作ることを表す．
 */
struct Cross {
    int x, y;
    char S[52];
};

/**
 * 答えがYESのときyn = 1, NOのときyn = 0で呼び出す．
 * mは交配の回数を表す．
 * cは交配の情報を持ったcross型の配列．
 */
void Output(int yn, int m, struct Cross* c) {
    if (!yn) {
        printf("NO\n");
    }
    else {
        printf("YES\n");
        printf("%d\n", m);
        int i;
        for (i = 0; i < m; i++)
            printf("%d %d %s\n", c[i].x, c[i].y, c[i].S);
    }
    fflush(stdout);
}


constexpr int MAX_N = 2005;
constexpr char YES[] = "YES";
constexpr char NO[] = "NO";
constexpr int UP = 0;
constexpr int LOW = 1;

class Selector {
public:
    Selector(const int n, const vector<string>& S, const vector<int> seeds[26][2])
        : n(n), S(S) {
        for (int i = 0; i < n; ++i) {
            isErasable[i] = true;
        }

        for (int i = 0; i < 26; ++i) {
            for (int j = 0; j < 2; ++j) {
                count[i][j] = seeds[i][j].size();
            }
        }

        for (int i = 0; i < 26; ++i) {
            for (int j = 0; j < 2; ++j) {
                this->seeds[i][j] = &seeds[i][j];
            }
        }
    }

    void Select(vector<int> result[26][2]) {
        for (int i = 0; i < 26; ++i) {
            RemoveIf(i);
        }

        for (int i = 0; i < 26; ++i) {
            for (int j = 0; j < 2; ++j) {
                result[i][j].reserve(count[i][j]);
                for (auto k : *seeds[i][j]) {
                    if (!isErasable[k]) continue;
                    result[i][j].push_back(k);
                }
            }
        }
    }

private:
    const int n;
    const vector<string>& S;
    const vector<int>* seeds[26][2];
    bool isErasable[MAX_N];
    int count[26][2];

    inline void RemoveIf(int index) {
        if ((count[index][UP] == 0) ^
            (count[index][LOW] == 0)) {
            const auto current = seeds[index];
            const auto target = count[index][UP] ? *current[UP] : *current[LOW];
            count[index][UP] = count[index][LOW] = 0;

            for (auto i : target) {
                if (!isErasable[i]) continue;
                isErasable[i] = false;

                for (auto c : S[i]) {
                    const int nx = toupper(c) - 'A';
                    if (nx == index) continue;
                    const int up = isupper(c) ? 0 : 1;
                    if (count[nx][up] == 0) continue;
                    count[nx][up]--;
                    assert(0 <= count[nx][up]);
                    RemoveIf(nx);
                }
            }
        }
    }
};


int n;
vector<string> S;
vector<int> seeds[26][2];

bool IsNo() {
    Selector selector(n, S, seeds);
    vector<int> result[26][2];
    selector.Select(result);

    for (int i = 0; i < 26; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (result[i][j].empty()) continue;
            return false;
        }
    }

    return true;
}

int main() {
    scanf("%d", &n);
    for (int i = 0; i < n; ++i) {
        string tmp;
        cin >> tmp;
        S.push_back(tmp);

        for (auto c : tmp) {
            const int isLow = islower(c) ? 1 : 0;
            const int index = toupper(c) - 'A';
            seeds[index][isLow].push_back(i);
        }
    }

    cout << (IsNo() ? NO : YES) << endl;
}
