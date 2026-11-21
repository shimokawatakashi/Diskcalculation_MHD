// eta_O, eta_H, eta_A [cm^2/s] table on a log-spaced (n_H, T, B) grid.
// usage (in Marchand2022/): make table, or ./eta_table ./input/input
// output columns: n_g T B eta_O eta_H eta_A  (n_g outer, B inner); n_g = gas particle density
// format = fortran: ii jj kk log10(n_H) log10(T) log10(B) eta_A eta_H eta_O as (3I4,6E13.5), no header

#include "typedef.hpp"
#include "calc_resistivity.hpp"
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::string Trim(const std::string &s)
{
    size_t b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return "";
    size_t e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

class TableInput
{
public:
    explicit TableInput(const std::string &file_name)
    {
        std::ifstream file(file_name);
        if (!file) {
            std::cerr << "error open file : " << file_name << std::endl;
            exit(1);
        }
        std::string line;
        while (std::getline(file, line)) {
            size_t c = line.find("//");
            if (c != std::string::npos) line = line.substr(0, c);
            line = Trim(line);
            if (line.empty() || line[0] == '#') continue;
            size_t eq = line.find('=');
            if (eq == std::string::npos) {
                std::cerr << "input: line without '=' : " << line << std::endl;
                exit(1);
            }
            kv_[Trim(line.substr(0, eq))] = Trim(line.substr(eq + 1));
        }
    }

    bool Has(const std::string &key) const { return kv_.count(key) > 0; }

    std::string Str(const std::string &key) const
    {
        if (!Has(key)) {
            std::cerr << "input: missing key '" << key << "'" << std::endl;
            exit(1);
        }
        return kv_.at(key);
    }
    std::string Str(const std::string &key, const std::string &def) const { return Has(key) ? kv_.at(key) : def; }
    double Num(const std::string &key) const { return std::stod(Str(key)); }
    double Num(const std::string &key, double def) const { return Has(key) ? std::stod(kv_.at(key)) : def; }
    int Int(const std::string &key) const { return std::stoi(Str(key)); }
    int Int(const std::string &key, int def) const { return Has(key) ? std::stoi(kv_.at(key)) : def; }

private:
    std::map<std::string, std::string> kv_;
};

DoubleVector LogSpace(const std::string &name, double vmin, double vmax, int num)
{
    if (!(vmin > 0.0) || !(vmax >= vmin) || num < 1 || (num == 1 && vmax != vmin)) {
        std::cerr << "input: invalid " << name << " range (min=" << vmin << ", max=" << vmax
                  << ", num=" << num << "); need 0 < min <= max, num >= 1 (num = 1 only if min = max)"
                  << std::endl;
        exit(1);
    }
    DoubleVector v(num);
    for (int k = 0; k < num; k++)
        v[k] = (num == 1) ? vmin : vmin * std::pow(vmax / vmin, double(k) / double(num - 1));
    return v;
}

// integral_{r0}^{r0*z} a^p da
double PowInt(double r0, double z, double p)
{
    if (std::fabs(p + 1.0) < 1.0e-12) return std::log(z);
    return std::pow(r0, p + 1.0) * (std::pow(z, p + 1.0) - 1.0) / (p + 1.0);
}

// dust per unit n_g: size <a>, mass from <a^3>, n_j/n_g (dust mass = mu m_p n_g d)
void SetupDust(const TableInput &in, DoubleVector &a, DoubleVector &m, DoubleVector &nd_per_ng)
{
    const std::string model = in.Str("dust_model");
    const double d   = in.Num("dust_to_gas", 1.0e-2);
    const double rho = in.Num("rho_grain", 2.9);
    const double mu  = in.Num("mu_gas", 2.3);
    const double rho_d_per_ng = mu * PROTON_MASS * d;

    if (model == "single") {
        const double as = in.Num("a_single");
        a.assign(1, as);
        m.assign(1, 4.0/3.0 * M_PI * rho * CUB(as));
        nd_per_ng.assign(1, rho_d_per_ng / m[0]);
    } else if (model == "mrn") {
        const double amin = in.Num("a_min"), amax = in.Num("a_max");
        const int    nb   = in.Int("n_bin");
        const double s    = in.Num("slope", -3.5);
        if (!(amin > 0.0) || !(amax > amin) || nb < 1) {
            std::cerr << "input: invalid MRN parameters (need 0 < a_min < a_max, n_bin >= 1)" << std::endl;
            exit(1);
        }
        const double z = std::pow(amax / amin, 1.0 / double(nb));
        const double mass_tot = PowInt(amin, amax / amin, s + 3.0);
        a.resize(nb); m.resize(nb); nd_per_ng.resize(nb);
        for (int j = 0; j < nb; j++) {
            const double r0 = amin * std::pow(z, double(j));
            const double I0 = PowInt(r0, z, s), I1 = PowInt(r0, z, s + 1.0), I3 = PowInt(r0, z, s + 3.0);
            a[j] = I1 / I0;
            m[j] = 4.0/3.0 * M_PI * rho * (I3 / I0);
            nd_per_ng[j] = rho_d_per_ng * (I3 / mass_tot) / m[j];
        }
    } else {
        std::cerr << "input: dust_model must be 'mrn' or 'single' (got '" << model << "')" << std::endl;
        exit(1);
    }
}

// Fortran E13.5: [-]0.dddddE+xx right-justified in 13 columns
std::string FortranE13_5(double v)
{
    char buf[32];
    if (v == 0.0) {
        std::snprintf(buf, sizeof(buf), "%13s", "0.00000E+00");
        return buf;
    }
    char sci[32];
    std::snprintf(sci, sizeof(sci), "%.4e", std::fabs(v)); // d.dddde+xx, rounded to 5 digits
    const std::string s(sci);
    const std::string digits = s.substr(0, 1) + s.substr(2, 4);
    const int ex = std::stoi(s.substr(s.find('e') + 1)) + 1;
    char body[32];
    if (std::abs(ex) <= 99)
        std::snprintf(body, sizeof(body), "%s0.%sE%c%02d", v < 0.0 ? "-" : "", digits.c_str(), ex < 0 ? '-' : '+', std::abs(ex));
    else
        std::snprintf(body, sizeof(body), "%s0.%s%c%03d", v < 0.0 ? "-" : "", digits.c_str(), ex < 0 ? '-' : '+', std::abs(ex));
    std::snprintf(buf, sizeof(buf), "%13s", body);
    return buf;
}

} // namespace

int main(int argc, char *argv[])
{
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <input file>" << std::endl;
        return 1;
    }
    const TableInput in(argv[1]);

    const DoubleVector ng = LogSpace("ng", in.Num("ng_min"), in.Num("ng_max"), in.Int("ng_num"));
    const DoubleVector Tg = LogSpace("T",  in.Num("T_min"),  in.Num("T_max"),  in.Int("T_num"));
    const DoubleVector Bg = LogSpace("B",  in.Num("B_min"),  in.Num("B_max"),  in.Int("B_num"));
    const double zeta = in.Num("zeta", 5.0e-17);
    const std::string outname = in.Str("output", "./table/eta_table.dat");
    const std::string format = in.Str("format", "plain");
    if (format != "plain" && format != "fortran") {
        std::cerr << "input: format must be 'plain' or 'fortran' (got '" << format << "')" << std::endl;
        return 1;
    }

    bool doH = false, doNa = false, doK = false;
    {
        std::stringstream ss(in.Str("species", "K"));
        std::string tok;
        while (std::getline(ss, tok, ',')) {
            tok = Trim(tok);
            if      (tok == "H")  doH  = true;
            else if (tok == "Na") doNa = true;
            else if (tok == "K")  doK  = true;
            else { std::cerr << "input: unknown species '" << tok << "' (use K, Na, H)" << std::endl; return 1; }
        }
        if (int(doH) + int(doNa) + int(doK) != 1) {
            std::cerr << "input: species must be exactly one of K, Na, H (got '" << in.Str("species", "K") << "')" << std::endl;
            return 1;
        }
    }

    DoubleVector a, m, nd_per_ng;
    SetupDust(in, a, m, nd_per_ng);

    const int nN = ng.size(), nT = Tg.size(), nB = Bg.size();
    std::vector<double> eta(3 * size_t(nN) * nT * nB, 0.0);
    auto idx = [&](int i, int j, int b) { return 3 * ((size_t(i) * nT + j) * nB + b); };

    auto solve = [&](Calc_Resistivity &cr, int i, int j) {
        DoubleVector nd(nd_per_ng.size());
        for (size_t k = 0; k < nd.size(); k++) nd[k] = nd_per_ng[k] * ng[i];
        cr.SetDebugContext(i, j, 0.0, 0.0, Tg[j], ng[i]);
        cr.SetGasParameter_For_Resistivity(zeta, Tg[j], ng[i]);
        cr.SetDustParameter_For_Resistivity(nd, m, a);
        cr.Calc_ChargeState();
        for (int b = 0; b < nB; b++) {
            cr.SetMagneticField_For_Resistivity(Bg[b]);
            cr.Calc_MagneticResistivity();
            size_t p = idx(i, j, b);
            eta[p] = cr.GetEtaO(); eta[p + 1] = cr.GetEtaH(); eta[p + 2] = cr.GetEtaA();
        }
    };

    // T_min row along n_g, then each n_g along T (previous solution as initial guess)
    std::vector<Calc_Resistivity> column(nN);
    {
        Calc_Resistivity cr;
        cr.SetActiveSpecies(doH, doNa, doK);
        for (int i = 0; i < nN; i++) {
            solve(cr, i, 0);
            column[i] = cr;
        }
    }
    std::cerr << "eta_table: T_min row done (" << nN << " densities)" << std::endl;

    for (int i = 0; i < nN; i++) {
        Calc_Resistivity cr = column[i];
        for (int j = 1; j < nT; j++) solve(cr, i, j);
    }

    const std::filesystem::path outdir = std::filesystem::path(outname).parent_path();
    if (!outdir.empty()) std::filesystem::create_directories(outdir);
    std::ofstream out(outname);
    if (!out) {
        std::cerr << "error open file : " << outname << std::endl;
        return 1;
    }
    if (format == "fortran") {
        // Fortran read(13,'(3I4,6E13.5)') ii, jj, kk, log10(n_g), log10(T), log10(B), etaA, etaH, etaO
        // (1-based indices, no header)
        char buf[32];
        for (int i = 0; i < nN; i++)
            for (int j = 0; j < nT; j++)
                for (int b = 0; b < nB; b++) {
                    size_t p = idx(i, j, b);
                    std::snprintf(buf, sizeof(buf), "%4d%4d%4d", i + 1, j + 1, b + 1);
                    out << buf << FortranE13_5(std::log10(ng[i])) << FortranE13_5(std::log10(Tg[j]))
                        << FortranE13_5(std::log10(Bg[b])) << FortranE13_5(eta[p + 2])
                        << FortranE13_5(eta[p + 1]) << FortranE13_5(eta[p]) << "\n";
                }
        std::cerr << "eta_table: wrote " << size_t(nN) * nT * nB << " rows (fortran 3I4,6E13.5) to "
                  << outname << std::endl;
        return 0;
    }
    out << "# eta table (Marchand et al. 2022 charge solver), input: " << argv[1] << "\n"
        << "# ng: " << ng.front() << " - " << ng.back() << " (" << nN << " log points) [cm^-3, gas particles; n_H = n_g/0.585]\n"
        << "# T : " << Tg.front() << " - " << Tg.back() << " (" << nT << " log points) [K]\n"
        << "# B : " << Bg.front() << " - " << Bg.back() << " (" << nB << " log points) [G]\n"
        << "# zeta = " << zeta << " s^-1, species = " << in.Str("species", "K")
        << ", dust_model = " << in.Str("dust_model") << ", dust bins = " << a.size()
        << ", grains: charge-resolved, Epstein drag (Tsukamoto & Okuzumi 2022)\n"
        << "# columns: ng[cm^-3] T[K] B[G] eta_O[cm^2/s] eta_H[cm^2/s] eta_A[cm^2/s]"
        << "  (order: ng outer, T, B inner)\n";
    out << std::scientific << std::setprecision(6);
    for (int i = 0; i < nN; i++)
        for (int j = 0; j < nT; j++)
            for (int b = 0; b < nB; b++) {
                size_t p = idx(i, j, b);
                out << ng[i] << " " << Tg[j] << " " << Bg[b] << " "
                    << eta[p] << " " << eta[p + 1] << " " << eta[p + 2] << "\n";
            }
    std::cerr << "eta_table: wrote " << size_t(nN) * nT * nB << " rows to " << outname << std::endl;
    return 0;
}
