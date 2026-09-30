// For WINDOWS
// Battery_sort 1.1.0
// Нове/New: вибір 0/1 для внутрішнього опору; відбір комірок з мінімальним розкидом;
//           розподіл по групах з вирівнюванням сум ємності (та опору) між групами.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
#include <windows.h>

using namespace std;

struct Cell { int id; double c; double r; };   // c - mAh, r - mOhm

// Вага опору відносно ємності (ємність - головний критерій, опір - другорядний).
// 0 = ігнорувати опір, 1 = опір рівноцінний ємності.
static const double W_IR = 0.1;

static double readNum(const string& prompt, double minV, bool integer) {
    while (true) {
        cout << prompt;
        string s;
        if (!getline(cin, s)) exit(1);
        replace(s.begin(), s.end(), ',', '.');
        try {
            size_t pos = 0;
            double v = stod(s, &pos);
            while (pos < s.size() && isspace((unsigned char)s[pos])) pos++;
            if (pos != s.size()) throw 1;
            if (v < minV || (integer && v != floor(v))) throw 1;
            return v;
        } catch (...) {
            cout << "Некоректне значення, повторіть.\n";
        }
    }
}

static void waitExit() {
    system("pause");
}

// ---------- відбір M комірок з мінімальним розкидом ----------
struct Stat { double sc = 0, scc = 0, sr = 0, srr = 0; int n = 0; };

static double spreadCost(const Stat& t, bool useR) {
    double n = t.n, mc = t.sc / n;
    double j = sqrt(max(0.0, t.scc / n - mc * mc)) / mc;
    if (useR) {
        double mr = t.sr / n;
        j += W_IR * sqrt(max(0.0, t.srr / n - mr * mr)) / mr;
    }
    return j;
}

static void addCell(Stat& t, const Cell& c, int sign) {
    t.sc += sign * c.c;  t.scc += sign * c.c * c.c;
    t.sr += sign * c.r;  t.srr += sign * c.r * c.r;
    t.n += sign;
}

static vector<Cell> selectCells(vector<Cell> all, int M, bool useR) {
    int N = all.size();
    if (M == N) return all;
    sort(all.begin(), all.end(), [](const Cell& a, const Cell& b) { return a.c < b.c; });
    // Мультистарт: кожне вікно з M сусідніх за ємністю (без опору це і є оптимум),
    // з опором - додатково обміни вибрана <-> невибрана, доки зменшується розкид.
    vector<char> bestSel; double bestJ = 1e300;
    for (int st = 0; st + M <= N; st++) {
        vector<char> sel(N, 0);
        Stat t;
        for (int k = st; k < st + M; k++) { sel[k] = 1; addCell(t, all[k], +1); }
        if (useR) {
            for (int it = 0; it < 100000; it++) {
                double best = spreadCost(t, true) - 1e-12;
                int bi = -1, bj = -1;
                for (int i = 0; i < N; i++) if (sel[i])
                    for (int j = 0; j < N; j++) if (!sel[j]) {
                        Stat u = t;
                        addCell(u, all[i], -1); addCell(u, all[j], +1);
                        double jc = spreadCost(u, true);
                        if (jc < best) { best = jc; bi = i; bj = j; }
                    }
                if (bi < 0) break;
                addCell(t, all[bi], -1); addCell(t, all[bj], +1);
                sel[bi] = 0; sel[bj] = 1;
            }
        }
        double jt = spreadCost(t, useR);
        if (jt < bestJ - 1e-12) { bestJ = jt; bestSel = sel; }
    }
    vector<char> sel = bestSel;
    vector<Cell> out;
    for (int i = 0; i < N; i++) if (sel[i]) out.push_back(all[i]);
    return out;
}

// ---------- розподіл по групах ----------
static double groupCost(const vector<double>& C, const vector<double>& Y, bool useR) {
    int s = C.size();
    double mn = 1e300, mx = -1e300, sum = 0, sq = 0;
    for (double c : C) { mn = min(mn, c); mx = max(mx, c); sum += c; sq += c * c; }
    double m = sum / s;
    double j = (mx - mn) / m + sqrt(max(0.0, sq / s - m * m)) / m;
    if (useR) {
        mn = 1e300; mx = -1e300; sum = 0; sq = 0;
        for (double y : Y) { double r = 1.0 / y; mn = min(mn, r); mx = max(mx, r); sum += r; sq += r * r; }
        m = sum / s;
        j += W_IR * ((mx - mn) / m + sqrt(max(0.0, sq / s - m * m)) / m);
    }
    return j;
}

static vector<vector<Cell>> buildGroups(vector<Cell> cells, int s, int p, bool useR) {
    sort(cells.begin(), cells.end(), [](const Cell& a, const Cell& b) { return a.c > b.c; });
    vector<vector<Cell>> g(s);
    for (size_t k = 0; k < cells.size(); k++) {          // "змійка"
        int row = k / s, col = k % s;
        g[(row % 2 == 0) ? col : s - 1 - col].push_back(cells[k]);
    }
    vector<double> C(s, 0), Y(s, 0);
    for (int i = 0; i < s; i++) for (auto& c : g[i]) { C[i] += c.c; if (useR) Y[i] += 1.0 / c.r; }
    double cur = groupCost(C, Y, useR);
    for (int pass = 0; pass < 2000; pass++) {
        bool improved = false;
        for (int a = 0; a < s; a++) for (int b = a + 1; b < s; b++)
            for (int i = 0; i < p; i++) for (int j = 0; j < p; j++) {
                Cell &x = g[a][i], &y = g[b][j];
                double oC1 = C[a], oC2 = C[b], oY1 = Y[a], oY2 = Y[b];
                C[a] += y.c - x.c; C[b] += x.c - y.c;
                if (useR) { Y[a] += 1.0 / y.r - 1.0 / x.r; Y[b] += 1.0 / x.r - 1.0 / y.r; }
                double nc = groupCost(C, Y, useR);
                if (nc < cur - 1e-12) { cur = nc; swap(x, y); improved = true; }
                else { C[a] = oC1; C[b] = oC2; Y[a] = oY1; Y[b] = oY2; }
            }
        if (!improved) break;
    }
    return g;
}

int main() {
    SetConsoleOutputCP(65001);
    SetConsoleCP(65001);
    cout << "Програма сортування та групування акумуляторів (Li-ion / LiFePO4)\n\n";

    cout << "Виберіть тип акумуляторів:\n";
    cout << "1 - Li-ion (мін 3.0 V, ном 3.7 V, макс 4.2 V)\n";
    cout << "2 - LiFePO4 (мін 2.5 V, ном 3.2 V, макс 3.65 V)\n";
    int chem_type = (int)readNum("Ваш вибір (1 або 2): ", 1, true);
    while (chem_type != 1 && chem_type != 2) chem_type = (int)readNum("Введіть 1 або 2: ", 1, true);

    double v_min_per = (chem_type == 1 ? 3.0 : 2.5);
    double v_nom_per = (chem_type == 1 ? 3.7 : 3.2);
    double v_max_per = (chem_type == 1 ? 4.2 : 3.65);
    string chem_name = (chem_type == 1 ? "Li-ion" : "LiFePO4");

    int N_total = (int)readNum("\nВведіть загальну кількість акумуляторів (всі доступні): ", 1, true);
    int M = (int)readNum("Скільки з них використати для збірки (рекомендується менше або дорівнює загальній): ", 1, true);
    int s = (int)readNum("Введіть розмір (довжина s — кількість секцій послідовно): ", 1, true);
    int p = (int)readNum("Введіть розмір (висота p — кількість паралельних комірок в групі): ", 1, true);

    if (M != s * p || N_total < M) {
        cout << "\nПОМИЛКА: Кількість комірок для збірки (M) повинна дорівнювати s * p.\n";
        cout << "Також загальна кількість повинна бути не меншою за M.\n";
        waitExit();
        return 1;
    }

    bool useR = readNum("\nВводити внутрішній опір комірок? 1 - так, 0 - ні: ", 0, true) == 1;

    cout << "\nВведення ємності акумуляторів (в mAh)" << (useR ? " та внутрішнього опору (в mΩ)" : "") << "\n";
    vector<Cell> all;
    for (int i = 1; i <= N_total; ++i) {
        Cell c; c.id = i; c.r = 0;
        c.c = readNum("Ємність акумулятора #" + to_string(i) + " (mAh): ", 1e-9, false);
        if (useR) c.r = readNum("Внутрішній опір акумулятора #" + to_string(i) + " (mΩ): ", 1e-9, false);
        all.push_back(c);
    }

    vector<Cell> chosen = selectCells(all, M, useR);
    vector<vector<Cell>> groups = buildGroups(chosen, s, p, useR);

    // ---------- звіт ----------
    ostringstream o;
    o << fixed;
    o << "ГРУПИ АКУМУЛЯТОРІВ (для паралельного з'єднання, " << chem_name << ")\n";
    o << "Вибрано " << M << " акумуляторів з " << N_total << " з мінімальним розкидом "
      << (useR ? "(ємність + опір)" : "(ємність)") << ".\n";
    o << (useR ? "Групи підібрані так, щоб суми ємності та опір груп були максимально близькими." : "Групи підібрані так, щоб суми ємності груп були максимально близькими.") << "\n";
    vector<char> used(N_total + 1, 0);
    for (auto& c : chosen) used[c.id] = 1;
    o << "Не використано: ";
    bool any = false;
    for (auto& c : all) if (!used[c.id]) { o << (any ? ", " : "") << "#" << c.id; any = true; }
    if (!any) o << "-";
    o << "\n\n";

    vector<double> gc(s), gr(s);
    for (int i = 0; i < s; i++) {
        double sum = 0, y = 0;
        o << "Група " << i + 1 << ": ";
        sort(groups[i].begin(), groups[i].end(), [](const Cell& a, const Cell& b) { return a.id < b.id; });
        for (auto& c : groups[i]) {
            sum += c.c; if (useR) y += 1.0 / c.r;
            o << "#" << c.id << " (" << setprecision(0) << c.c << " mAh";
            if (useR) o << ", " << setprecision(1) << c.r << " mΩ";
            o << ") ";
        }
        gc[i] = sum; gr[i] = useR ? 1.0 / y : 0;
        o << "| сер. = " << setprecision(1) << sum / p << " mAh | сума = " << setprecision(0) << sum << " mAh";
        if (useR) o << " | R = " << setprecision(2) << gr[i] << " mΩ";
        o << "\n";
    }
    double cmin = *min_element(gc.begin(), gc.end()), cmax = *max_element(gc.begin(), gc.end());
    double cavg = 0; for (double c : gc) cavg += c; cavg /= s;
    o << "\nРОЗКИД МІЖ ГРУПАМИ\n";
    o << "Ємність: " << setprecision(0) << cmax - cmin << " mAh (" << setprecision(2) << (cmax - cmin) / cavg * 100 << " %)\n";
    if (useR) {
        double rmin = *min_element(gr.begin(), gr.end()), rmax = *max_element(gr.begin(), gr.end());
        double rsum = 0; for (double r : gr) rsum += r;
        o << "Опір груп: " << setprecision(2) << rmin << " ... " << rmax << " mΩ ("
          << (rmax - rmin) / (rsum / s) * 100 << " %)\n";
        o << "Опір збірки (без шин): " << rsum << " mΩ\n";
    }

    double min_v = v_min_per * s, nom_v = v_nom_per * s, max_v = v_max_per * s;
    o << "\nНАПРУГА ЗБІРКИ (" << s << "S, " << chem_name << ")\n" << setprecision(2);
    o << "Мінімальна:   " << min_v << " V\n";
    o << "Робоча (номінал): " << nom_v << " V\n";
    o << "Максимальна:  " << max_v << " V\n";

    double pack_ah = cmin / 1000.0;           // ліміт - найслабша група
    double pack_wh = nom_v * pack_ah;
    o << "\nПРИБЛИЗНА ЄМНІСТЬ ЗБІРКИ\n";
    o << "Ємність в ампер-годинах:   " << setprecision(3) << pack_ah << " Ah (по найслабшій групі)\n";
    o << "Ємність в ват-годинах:    " << setprecision(1) << pack_wh << " Wh\n";
    o << "(Середня ємність групи " << setprecision(0) << cavg << " mAh)\n";

    cout << "\n" << o.str();

    int save_choice = (int)readNum("\nЗберегти результати обчислень у текстовий файл? 1 - так , 0 - ні: ", 0, true);
    if (save_choice == 1) {
        time_t t_now = time(nullptr);
        char b1[64], b2[64];
        tm* tm_local = localtime(&t_now);
        strftime(b1, sizeof b1, "%Y-%m-%d_%H-%M-%S", tm_local);
        strftime(b2, sizeof b2, "%Y-%m-%d %H:%M:%S", tm_local);
        string filename = string("Battery_Sort_") + b1 + ".txt";
        ofstream fout(filename);
        if (fout.is_open()) {
            fout << "РЕЗУЛЬТАТИ СОРТУВАННЯ ТА ГРУПУВАННЯ АКУМУЛЯТОРІВ (" << chem_name << ")\n\n";
            fout << "Загальна кількість введених акумуляторів: " << N_total << "\n";
            fout << "Використано для збірки: " << M << "\n";
            fout << "Конфігурація: " << s << "S " << p << "P\n\n";
            fout << o.str();
            fout << "\nФайл згенеровано автоматично " << b2 << "\n";
            fout.close();
            cout << "\nРезультати успішно збережено у файл: " << filename << "\n";
        } else {
            cout << "\nНе вдалося створити файл!\n";
        }
    }

    cout << "\nПрограма завершена.\n";
    waitExit();
    return 0;
}
