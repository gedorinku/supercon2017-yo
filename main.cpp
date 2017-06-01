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
void Output(int yn, int m, struct Cross *c) {
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


constexpr int MAX_N = 2005;
constexpr char YES[] = "YES";
constexpr char NO[] = "NO";
constexpr int UP = 0;
constexpr int LOW = 1;

int n;
vector<string> S;
unordered_set<int> seeds[26][2];

inline void RemoveIf(int index) {
    auto &current = seeds[index];
    if (seeds[index][UP].empty() ^
        seeds[index][LOW].empty()) {
        auto &target = current[0].empty() ? current[1] : current[0];
        for (auto i : target) {
            target.erase(i);
            for (auto c : S[i]) {
                int nx = toupper(c) - 'A';
                if (nx == index) continue;
                RemoveIf(nx);
            }
            if (target.empty()) break;
        }
   }
}

bool IsNo() {
    for (int i = 0; i < 26; ++i) {
        RemoveIf(i);
    }

    for (int i = 0; i < 26; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (seeds[i][j].empty()) continue;
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
            cout << index << " " << isLow << endl;
            seeds[index][isLow].insert(i);
        }
    }

    cout << (IsNo() ? NO : YES) << endl;
}
