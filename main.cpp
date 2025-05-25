#include <iostream>

#include <string>

#include <vector>

#include <memory>

#include <algorithm>

#include <cctype>

#include <stdexcept>

#include <map>

#include <iomanip>

#include <sstream>

#include <cmath>

#include <set>

using namespace std;

inline string to_lower_util(string s) {
    transform(s.begin(), s.end(), s.begin(),
              [](unsigned char c) {
                  return tolower(c);
              });
    return s;
}

string trim_string_util(const string & str) {
    const string whitespace = " \t\n\r\f\v";
    size_t start = str.find_first_not_of(whitespace);
    if (start == string::npos) return "";
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

double parse_value_with_metric_prefix_util(const string & val_str_orig) {
    string val_str = trim_string_util(val_str_orig);
    if (val_str.empty()) throw invalid_argument("Empty value string for parsing.");
    size_t suffix_start_pos = 0;
    double base_val;
    try {
        base_val = stod(val_str, & suffix_start_pos);
    } catch (const exception & ) {
        throw invalid_argument("Invalid numeric value in string: \"" + val_str_orig + "\"");
    }
    if (suffix_start_pos >= val_str.length()) return base_val;
    string suffix_str = to_lower_util(val_str.substr(suffix_start_pos));
    double multiplier = 1.0;
    if (suffix_str.rfind("meg", 0) == 0) multiplier = 1e6;
    else if (!suffix_str.empty()) {
        char prefix_char = suffix_str[0];
        switch (prefix_char) {
            case 'p':
                multiplier = 1e-12;
                break;
            case 'n':
                multiplier = 1e-9;
                break;
            case 'u':
                multiplier = 1e-6;
                break;
            case 'm':
                multiplier = 1e-3;
                break;
            case 'k':
                multiplier = 1e3;
                break;
            case 'g':
                multiplier = 1e9;
                break;
        }
    }
    return base_val * multiplier;
}

vector < double > gaussian_elimination_matrix(vector < vector < double >> A, vector < double > b) {
    int n = A.size();
    if (n == 0 || (n > 0 && (A[0].size() != n || b.size() != n))) throw runtime_error("Invalid matrix or vector dimensions for Gaussian elimination.");
    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k)
            if (abs(A[k][i]) > abs(A[max_row][i])) max_row = k;
        swap(A[i], A[max_row]);
        swap(b[i], b[max_row]);
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[i][i]) < 1e-12) continue;
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    vector < double > x(n);
    for (int i = n - 1; i >= 0; --i) {
        if (abs(A[i][i]) < 1e-12) {
            if (abs(b[i]) > 1e-12) throw runtime_error("System is inconsistent (no solution).");
            x[i] = 0;
        } else {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j) x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
    }
    return x;
}

class Circuit;

class Component {
public: string name;
    string node1_name,
            node2_name;
    double value;
    Component(string name_val, string n1, string n2, double val = 0.0): name(move(name_val)),
                                                                        node1_name(move(n1)),
                                                                        node2_name(move(n2)),
                                                                        value(val) {}
    Component(string name_val, string n1, string n2,
              const string & val_str): name(move(name_val)),
                                       node1_name(move(n1)),
                                       node2_name(move(n2)) {
        this -> value = parse_value_with_metric_prefix_util(val_str);
    }
    virtual~Component() =
    default;
    virtual void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m_map, double h,
                       const vector < double > & prev_sol) = 0;
    virtual string to_netlist_string() const = 0;
    virtual void update_time_dependant_value(double time) {}
};

using ResultPoint = map < string, double > ;

class VoltageSource;
class Inductor;

class Circuit {
public: vector < unique_ptr < Component >> components;
    map < string,
            int > node_to_idx;
    vector < string > idx_to_node_name;
    vector < VoltageSource * > voltage_source_list;
    vector < Inductor * > inductor_list;
    string ground_node_explicit_name = "0";
    bool ground_node_exists = true;
    bool tran_solved = false;
    vector < ResultPoint > tran_results;

    void add_component(unique_ptr < Component > comp);
    bool is_ground(const string & node_name) const {
        return node_name == ground_node_explicit_name;
    }
    int get_node_matrix_index(const string & node_name) const;
    int prepare_for_analysis();
    void build_mna_matrix(vector < vector < double >> & A, vector < double > & z, double h,
                          const vector < double > & prev_sol);
    void perform_transient_analysis(double t_step, double t_stop);
    void clear();
};

class Resistor: public Component {
public: Resistor(const string & name,
                 const string & n1,
                 const string & n2,
                 const string & val): Component(name, n1, n2, val) {}
    string to_netlist_string() const override;
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
               const vector < double > & p) override;
};

class Capacitor: public Component {
public: Capacitor(const string & name,
                  const string & n1,
                  const string & n2,
                  const string & val): Component(name, n1, n2, val) {}
    string to_netlist_string() const override;
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
               const vector < double > & p) override;
};

class Inductor: public Component {
public: Inductor(const string & name,
                 const string & n1,
                 const string & n2,
                 const string & val): Component(name, n1, n2, val) {}
    string to_netlist_string() const override;
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
               const vector < double > & p) override;
};

class VoltageSource: public Component {
public: enum SourceType {
        DC,
        SINUSOIDAL,
        PULSE
    };
    SourceType sourceType;
    double dc_offset,
            amplitude,
            frequency;
    double v1,
            v2,
            td,
            tr,
            tf,
            pw,
            per;
    vector < string > raw_params;

    VoltageSource(const string & v_name,
                  const string & n1,
                  const string & n2,
                  const vector < string > & params);
    string to_netlist_string() const override;
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
               const vector < double > & p) override;
    void update_time_dependant_value(double time) override;
};

class CurrentSource: public Component {
public: CurrentSource(const string & name,
                      const string & n1,
                      const string & n2,
                      const string & val): Component(name, n1, n2, val) {}
    string to_netlist_string() const override;
    void stamp(Circuit & circuit, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
               const vector < double > & p) override;
};

void Circuit::clear() {
    components.clear();
    node_to_idx.clear();
    idx_to_node_name.clear();
    voltage_source_list.clear();
    inductor_list.clear();
    ground_node_exists = true;
    ground_node_explicit_name = "0";
    tran_solved = false;
    tran_results.clear();
}

void Circuit::add_component(unique_ptr < Component > comp) {
    components.push_back(move(comp));
    tran_solved = false;
}

int Circuit::get_node_matrix_index(const string & node_name) const {
    if (is_ground(node_name)) return -1;
    auto it = node_to_idx.find(node_name);
    return (it != node_to_idx.end()) ? it -> second : -2;
}

int Circuit::prepare_for_analysis() {
    node_to_idx.clear();
    idx_to_node_name.clear();
    voltage_source_list.clear();
    inductor_list.clear();
    set < string > unique_node_names;
    for (const auto & comp: components) {
        unique_node_names.insert(comp -> node1_name);
        unique_node_names.insert(comp -> node2_name);
        if (auto vs = dynamic_cast < VoltageSource * > (comp.get())) voltage_source_list.push_back(vs);
        else if (auto ind = dynamic_cast < Inductor * > (comp.get())) inductor_list.push_back(ind);
    }
    int idx = 0;
    for (const auto & name: unique_node_names) {
        if (!is_ground(name)) {
            node_to_idx[name] = idx++;
            idx_to_node_name.push_back(name);
        }
    }
    return idx;
}

void Circuit::build_mna_matrix(vector < vector < double >> & A, vector < double > & z, double h,
                               const vector < double > & prev_sol) {
    if (!components.empty() && !ground_node_exists) throw runtime_error("No ground node defined.");
    int N = prepare_for_analysis();
    int M = voltage_source_list.size() + inductor_list.size();
    int system_size = N + M;
    if (system_size == 0) {
        A.clear();
        z.clear();
        return;
    }

    vector < vector < double >> G(N, vector < double > (N, 0.0)), B(N, vector < double > (M, 0.0));
    vector < vector < double >> C(M, vector < double > (N, 0.0)), D(M, vector < double > (M, 0.0));
    vector < double > J(N, 0.0), E(M, 0.0);
    map < string, int > m_map;
    int m_counter = 0;
    for (const auto & vs: voltage_source_list) m_map[vs -> name] = m_counter++;
    for (const auto & l: inductor_list) m_map[l -> name] = m_counter++;

    for (const auto & comp: components) {
        comp -> stamp( * this, G, B, C, D, J, E, m_map, h, prev_sol);
    }

    A.assign(system_size, vector < double > (system_size, 0.0));
    z.assign(system_size, 0.0);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) A[i][j] = G[i][j];
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) A[i][N + j] = B[i][j];
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j) A[N + i][j] = C[i][j];
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) A[N + i][N + j] = D[i][j];
    for (int i = 0; i < N; ++i) z[i] = J[i];
    for (int i = 0; i < M; ++i) z[N + i] = E[i];
}

void Circuit::perform_transient_analysis(double t_step, double t_stop) {
    tran_solved = false;
    tran_results.clear();
    if (t_step <= 0 || t_stop <= 0 || t_step > t_stop) throw runtime_error("Invalid transient parameters.");

    int N = prepare_for_analysis();
    int M = voltage_source_list.size() + inductor_list.size();
    if (N + M == 0) {
        tran_solved = true;
        return;
    }

    vector < double > prev_solution(N + M, 0.0);

    for (double t = 0; t <= t_stop + (t_step / 2.0); t += t_step) {
        for (auto & comp: components) {
            comp -> update_time_dependant_value(t);
        }
        vector < vector < double >> A;
        vector < double > z;
        build_mna_matrix(A, z, t_step, prev_solution);

        vector < double > current_solution = gaussian_elimination_matrix(A, z);

        ResultPoint result_at_t;
        result_at_t["time"] = t;
        for (int i = 0; i < N; ++i) result_at_t["V(" + idx_to_node_name[i] + ")"] = current_solution[i];

        map < string, int > m_map;
        int m_counter = 0;
        for (const auto & vs: voltage_source_list) m_map[vs -> name] = m_counter++;
        for (const auto & l: inductor_list) m_map[l -> name] = m_counter++;
        for (const auto & p: m_map) result_at_t["I(" + p.first + ")"] = current_solution[N + p.second];

        tran_results.push_back(result_at_t);
        prev_solution = current_solution;
    }
    tran_solved = true;
    cout << "Transient analysis completed." << endl;
}

string Resistor::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}
void Resistor::stamp(Circuit & c, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
                     const vector < double > & p) {
    double g = 1.0 / value;
    int i1 = c.get_node_matrix_index(node1_name);
    int i2 = c.get_node_matrix_index(node2_name);
    if (i1 >= 0) G[i1][i1] += g;
    if (i2 >= 0) G[i2][i2] += g;
    if (i1 >= 0 && i2 >= 0) {
        G[i1][i2] -= g;
        G[i2][i1] -= g;
    }
}

string Capacitor::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}
void Capacitor::stamp(Circuit & c, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
                      const vector < double > & p) {
    if (h == 0) return;
    double g_eq = value / h;
    int i1 = c.get_node_matrix_index(node1_name);
    int i2 = c.get_node_matrix_index(node2_name);
    if (i1 >= 0) G[i1][i1] += g_eq;
    if (i2 >= 0) G[i2][i2] += g_eq;
    if (i1 >= 0 && i2 >= 0) {
        G[i1][i2] -= g_eq;
        G[i2][i1] -= g_eq;
    }
    if (!p.empty()) {
        double v1_prev = (i1 >= 0) ? p[i1] : 0.0;
        double v2_prev = (i2 >= 0) ? p[i2] : 0.0;
        double i_eq = g_eq * (v1_prev - v2_prev);
        if (i1 >= 0) J[i1] += i_eq;
        if (i2 >= 0) J[i2] -= i_eq;
    }
}

string Inductor::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}
void Inductor::stamp(Circuit & c, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
                     const vector < double > & p) {
    int m_idx = m.at(name);
    int i1 = c.get_node_matrix_index(node1_name);
    int i2 = c.get_node_matrix_index(node2_name);
    if (i1 >= 0) {
        B[i1][m_idx] += 1.0;
        C[m_idx][i1] += 1.0;
    }
    if (i2 >= 0) {
        B[i2][m_idx] -= 1.0;
        C[m_idx][i2] -= 1.0;
    }
    D[m_idx][m_idx] -= (h > 0 ? value / h : 1e-12);
    if (h > 0 && !p.empty()) {
        int N = c.idx_to_node_name.size();
        E[m_idx] -= (value / h) * p[N + m_idx];
    }
}

VoltageSource::VoltageSource(const string & v_name,
                             const string & n1,
                             const string & n2,
                             const vector < string > & params): Component(v_name, n1, n2, "0"), sourceType(DC), dc_offset(0), amplitude(0), frequency(0), v1(0), v2(0), td(0), tr(0), tf(0), pw(0), per(0), raw_params(params) {
    if (params.empty()) throw runtime_error("No value/params for voltage source " + name);
    string first_param_upper = to_lower_util(params[0]);
    if (first_param_upper.rfind("sin", 0) == 0) {
        string combined_params;
        for (size_t i = 0; i < params.size(); ++i) {
            combined_params += params[i] + " ";
        }
        size_t start_paren = combined_params.find('('), end_paren = combined_params.rfind(')');
        if (start_paren == string::npos || end_paren == string::npos) throw runtime_error("Mismatched parentheses in SIN() for " + name);
        string sin_params_str = combined_params.substr(start_paren + 1, end_paren - start_paren - 1);
        stringstream ss(sin_params_str);
        string offset_str, amp_str, freq_str;
        ss >> offset_str >> amp_str >> freq_str;
        if (ss.fail() || !ss.eof()) throw runtime_error("Invalid SIN() params for " + name);
        sourceType = SINUSOIDAL;
        dc_offset = parse_value_with_metric_prefix_util(offset_str);
        amplitude = parse_value_with_metric_prefix_util(amp_str);
        frequency = parse_value_with_metric_prefix_util(freq_str);
        this -> value = dc_offset;
    } else if (first_param_upper.rfind("pulse", 0) == 0) {
        string combined_params;
        for (size_t i = 0; i < params.size(); ++i) {
            combined_params += params[i] + " ";
        }
        size_t start_paren = combined_params.find('('), end_paren = combined_params.rfind(')');
        if (start_paren == string::npos || end_paren == string::npos) throw runtime_error("Mismatched parentheses in PULSE() for " + name);
        string pulse_params_str = combined_params.substr(start_paren + 1, end_paren - start_paren - 1);
        stringstream ss(pulse_params_str);
        string v1_s, v2_s, td_s, tr_s, tf_s, pw_s, per_s;
        ss >> v1_s >> v2_s >> td_s >> tr_s >> tf_s >> pw_s >> per_s;
        if (ss.fail() || !ss.eof()) throw runtime_error("Invalid PULSE() params for " + name);
        sourceType = PULSE;
        v1 = parse_value_with_metric_prefix_util(v1_s);
        v2 = parse_value_with_metric_prefix_util(v2_s);
        td = parse_value_with_metric_prefix_util(td_s);
        tr = parse_value_with_metric_prefix_util(tr_s);
        tf = parse_value_with_metric_prefix_util(tf_s);
        pw = parse_value_with_metric_prefix_util(pw_s);
        per = parse_value_with_metric_prefix_util(per_s);
        this -> value = v1;
    } else if (params.size() == 1) {
        sourceType = DC;
        this -> value = parse_value_with_metric_prefix_util(params[0]);
        this -> dc_offset = this -> value;
    } else throw runtime_error("Invalid params for voltage source " + name);
}

string VoltageSource::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " ";
    for (size_t i = 0; i < raw_params.size(); ++i) oss << raw_params[i] << (i == raw_params.size() - 1 ? "" : " ");
    return oss.str();
}

void VoltageSource::stamp(Circuit & c, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
                          const vector < double > & p) {
    int m_idx = m.at(name);
    int i1 = c.get_node_matrix_index(node1_name);
    int i2 = c.get_node_matrix_index(node2_name);
    if (i1 >= 0) {
        B[i1][m_idx] += 1.0;
        C[m_idx][i1] += 1.0;
    }
    if (i2 >= 0) {
        B[i2][m_idx] -= 1.0;
        C[m_idx][i2] -= 1.0;
    }
    E[m_idx] += this -> value;
}

void VoltageSource::update_time_dependant_value(double time) {
    if (sourceType == SINUSOIDAL) {
        if (frequency > 0) value = dc_offset + amplitude * sin(2 * 3.14159265358979323846 * frequency * time);
        else value = dc_offset + amplitude;
    } else if (sourceType == PULSE) {
        if (time < td) value = v1;
        else {
            double t_rel = per > 0 ? fmod(time - td, per) : time - td;
            if (tr > 0 && t_rel < tr) value = v1 + (v2 - v1) * (t_rel / tr);
            else if (t_rel < tr + pw) value = v2;
            else if (tf > 0 && t_rel < tr + pw + tf) value = v2 + (v1 - v2) * ((t_rel - tr - pw) / tf);
            else value = v1;
        }
    }
}

string CurrentSource::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}
void CurrentSource::stamp(Circuit & c, vector < vector < double >> & G, vector < vector < double >> & B, vector < vector < double >> & C, vector < vector < double >> & D, vector < double > & J, vector < double > & E, map < string, int > & m, double h,
                          const vector < double > & p) {
    int i1 = c.get_node_matrix_index(node1_name);
    int i2 = c.get_node_matrix_index(node2_name);
    if (i1 >= 0) J[i1] -= value;
    if (i2 >= 0) J[i2] += value;
}

struct Command {
    string type;
    vector < string > args;
};

class Parser {
public: Command parse_line(const string & line);
    void execute_command(Command & cmd, Circuit & circuit);
private: void handle_add_component(const vector < string > & args, Circuit & circuit);
    void handle_list_components(const Circuit & circuit);
    void handle_list_nodes(const Circuit & circuit);
};

Command Parser::parse_line(const string & line) {
    stringstream ss(trim_string_util(line));
    string segment;
    Command cmd;
    bool first = true;
    while (ss >> segment) {
        if (first) {
            cmd.type = segment;
            first = false;
        } else {
            cmd.args.push_back(segment);
        }
    }
    return cmd;
}

void Parser::execute_command(Command & cmd, Circuit & circuit) {
    try {
        string cmd_type_lower = to_lower_util(cmd.type);
        if (cmd_type_lower == "exit" || cmd_type_lower == "quit") {
            cout << "Exiting simulator." << endl;
            exit(0);
        } else if (cmd_type_lower == ".list") {
            handle_list_components(circuit);
        } else if (cmd_type_lower == ".nodes") {
            handle_list_nodes(circuit);
        } else if (cmd_type_lower == ".clear") {
            circuit.clear();
            cout << "Circuit cleared." << endl;
        } else {
            cmd.args.insert(cmd.args.begin(), cmd.type);
            handle_add_component(cmd.args, circuit);
        }
    } catch (const exception & e) {
        cerr << "Error: " << e.what() << endl;
    }
}

void Parser::handle_add_component(const vector < string > & args, Circuit & circuit) {
    if (args.size() < 4) throw runtime_error("Insufficient arguments. Usage: <name> <node1> <node2> <value> [params...]");
    string name = args[0];
    char type_char = toupper(name[0]);
    string n1 = args[1];
    string n2 = args[2];
    vector < string > remaining_args(args.begin() + 3, args.end());

    switch (type_char) {
        case 'R':
            circuit.add_component(make_unique < Resistor > (name, n1, n2, remaining_args[0]));
            break;
        case 'V':
            circuit.add_component(make_unique < VoltageSource > (name, n1, n2, remaining_args));
            break;
        case 'I':
            circuit.add_component(make_unique < CurrentSource > (name, n1, n2, remaining_args[0]));
            break;
        case 'C':
            circuit.add_component(make_unique < Capacitor > (name, n1, n2, remaining_args[0]));
            break;
        case 'L':
            circuit.add_component(make_unique < Inductor > (name, n1, n2, remaining_args[0]));
            break;
        default:
            throw runtime_error("Unknown component type: " + string(1, type_char));
    }
    cout << "Added component " << name << endl;
}

void Parser::handle_list_components(const Circuit & circuit) {
    cout << "Components in circuit:" << endl;
    if (circuit.components.empty()) {
        cout << "  (No components)" << endl;
        return;
    }
    for (const auto & comp: circuit.components) {
        cout << "  - " << comp -> to_netlist_string() << endl;
    }
}
void Parser::handle_list_nodes(const Circuit & circuit) {
    set < string > all_nodes;
    for (const auto & comp: circuit.components) {
        all_nodes.insert(comp -> node1_name);
        all_nodes.insert(comp -> node2_name);
    }
    cout << "Available nodes:" << endl;
    if (all_nodes.empty()) {
        cout << "  (No nodes in circuit)" << endl;
        return;
    }
    for (const auto & node_name: all_nodes) {
        cout << "  " << node_name << (circuit.is_ground(node_name) ? " (Ground)" : "") << endl;
    }
}

int main() {
    Circuit circuit_main;
    Parser parser_main;
    string line_input;
    cout << "Welcome to SUTSpice! Enter commands." << endl;
    while (true) {
        cout << "> ";
        if (!getline(cin, line_input)) break;
        if (line_input.empty()) continue;
        Command cmd = parser_main.parse_line(line_input);
        parser_main.execute_command(cmd, circuit_main);
    }
    return 0;
}