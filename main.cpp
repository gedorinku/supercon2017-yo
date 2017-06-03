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
constexpr int INF = 1 << 30;
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
                    const int up = isupper(c) ? UP : LOW;
                    if (count[nx][up] == 0) continue;
                    count[nx][up]--;
                    assert(0 <= count[nx][up]);
                    RemoveIf(nx);
                }
            }
        }
    }
};

struct Status {
    string seed;
    int eval;
    vector<int> process;
    vector<bool> used;
    vector<vector<bool>> usedAlpha;

    Status() : eval(INF) {
    }

    Status(int n) : eval(INF), used(n, false), usedAlpha(26, vector<bool>(2, false)) {
    }

    Status(const string& seed, const int eval, const vector<int>& process, const vector<bool>& used, const vector<vector<bool>>& usedAlpha)
        : seed(seed), eval(eval), process(process), used(used), usedAlpha(usedAlpha) {
    }

    bool operator<(const Status& right) const {
        return eval < right.eval;
    }

    bool operator>(const Status& right) const {
        return eval > right.eval;
    }
};


int n;
vector<string> S;
vector<int> seeds[26][2];
vector<int> selectedTable[26][2];
vector<int> selectedS;

bool IsNo() {
    Selector selector(n, S, seeds);
    selector.Select(selectedTable);

    for (int i = 0; i < 26; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (selectedTable[i][j].empty()) continue;
            for (auto k : selectedTable[i][j]) {
                selectedS.push_back(k);
            }
        }
    }

    if (selectedS.empty()) return true;

    sort(selectedS.begin(), selectedS.end());
    selectedS.erase(unique(selectedS.begin(), selectedS.end()), selectedS.end());

    return false;
}

inline int Evaluate(const string& seed) {
    int eval = 0;
    for (int i = 0; i < seed.size(); ++i) {
        if (isupper(seed[i])) {
            eval++;
            continue;
        }
        if (i == seed.size() - 1 || seed[i] != tolower(seed[i + 1])) {
            eval++;
            continue;
        }

        i++;
    }

    return eval;
}

inline void PushBackIfNotDuplicate(string& result, const char x) {
    if (result.empty() || (!result.empty() && result[result.size() - 1] != x)) {
        result.push_back(x);
    }
}

inline string Merge(const string& a, const string& b) {
    string result;

    for (int i = 0, j = 0; i < a.size() || j < b.size();) {
        if (a.size() <= i) {
            PushBackIfNotDuplicate(result, b[j]);
            j++;
            continue;
        }
        if (b.size() <= j) {
            PushBackIfNotDuplicate(result, a[i]);
            i++;
            continue;
        }

        const char ai = a[i], bj = b[j];
        const char lowA = tolower(ai), lowB = tolower(bj);
        if (lowA < lowB || (lowA == lowB && islower(ai))) {
            PushBackIfNotDuplicate(result, ai);
            i++;
        }
        else {
            PushBackIfNotDuplicate(result, bj);
            j++;
        }
    }

    return result;
}

void PrintAnswer(const Status& ans) {
    Cross cross[MAX_N];
    string current;
    {
        Cross c;
        c.x = ans.process[0] + 1;
        c.y = ans.process[1] + 1;
        current = Merge(S[c.x - 1], S[c.y - 1]);
        strcpy(c.S, current.c_str());
        cross[0] = c;
    }
    for (int i = 2; i < ans.process.size(); ++i) {
        Cross c;
        c.x = ans.process[i] + 1;
        c.y = n + i - 1;
        current = Merge(current, S[c.x - 1]);
        strcpy(c.S, current.c_str());
        cross[i - 1] = c;
    }
    {
        Cross c;
        c.x = c.y = n + ans.process.size() - 1;
        strcpy(c.S, "!");
        cross[ans.process.size() - 1] = c;
    }

    Output(1, ans.process.size(), cross);
}

void MarkAllKuse(vector<vector<bool>>& usedAlpha, const string& seed) {
    for (auto c : seed) {
        const int isLow = islower(c) ? LOW : UP;
        const int index = toupper(c) - 'A';
        usedAlpha[index][isLow] = true;
    }
}


Status current, best;

void Solve(int preUsed) {
    //cerr << "::" << preUsed << endl;
    vector<pair<int, int>> nextSeed;


    for (int i = 0; i < 26; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (current.usedAlpha[i][j]) continue;

            for (auto k : selectedTable[i][j]) {
                if (current.used[k] || k < preUsed) continue;

                nextSeed.push_back({ Evaluate(Merge(current.seed, S[k])), k });
            }
        }
    }

    sort(nextSeed.begin(), nextSeed.end());
    nextSeed.erase(unique(nextSeed.begin(), nextSeed.end()), nextSeed.end());
    //cerr << "::" << nextSeed.size() << endl;

    //for (auto p : nextSeed) {
    if (!nextSeed.empty()) {
        auto p = nextSeed[0];
        int i = p.second;
        string pre = current.seed;
        int preEval = current.eval;
        auto preUsedAlpha = current.usedAlpha;
        current.seed = Merge(current.seed, S[i]);
        current.eval = p.first;
        current.used[i] = true;
        //current.usedAlpha[i][j] = true;
        MarkAllKuse(current.usedAlpha, S[i]);

        current.process.push_back(i);

        bool finished = false;

        if (current.eval == 0) {
            if (best.process.empty() || current.process.size() < best.process.size()) {
                best = current;
                cerr << "best:" << best.process.size() << endl;
            }
            finished = true;
        }

        if (best.process.empty() || current.process.size() < best.process.size() - 1) {
            Solve(i);
        }

        current.seed = pre;
        current.eval = preEval;
        current.used[i] = false;
        current.usedAlpha = preUsedAlpha;
        current.process.pop_back();

        if (finished) return;
    }
}

int main() {
    auto start = chrono::high_resolution_clock::now();

    scanf("%d", &n);
    for (int i = 0; i < n; ++i) {
        string tmp;
        cin >> tmp;
        S.push_back(tmp);

        for (auto c : tmp) {
            const int isLow = islower(c) ? LOW : UP;
            const int index = toupper(c) - 'A';
            seeds[index][isLow].push_back(i);
        }
    }

    current = Status(n);
    best = Status(n);

    if (IsNo()) {
        Output(0, 0, nullptr);
    }
    else {
        vector<pair<int, int>> sortedS;
        for (auto i : selectedS) {
            sortedS.push_back(make_pair(Evaluate(S[i]), i));
        }
        sort(sortedS.begin(), sortedS.end());

        for (auto i : sortedS) {
            //cerr << i.second << "#" << endl;
            current = Status(n);

            current.seed = Merge(current.seed, S[i.second]);
            current.eval = Evaluate(current.seed);
            current.used[i.second] = true;
            MarkAllKuse(current.usedAlpha, current.seed);
            current.process.push_back(i.second);
            Solve(i.second);
        }
        PrintAnswer(best);
    }

    auto end = chrono::high_resolution_clock::now();
    auto dur = chrono::duration_cast<chrono::microseconds>(end - start);
    cerr << dur.count() / 1000.0 << "ms" << endl;

    int maxv = 0;
    for (int i = 0; i < 26; ++i) {
        for (int j = 0; j < 2; ++j) {
            if (selectedTable[i][j].empty()) continue;
            maxv++;
        }
    }
    cout << "max:" << maxv << endl;
}
