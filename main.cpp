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
              [](unsigned char c){ return tolower(c); });
    return s;
}

string trim_string_util(const string& str) {
    const string whitespace = " \t\n\r\f\v";
    size_t start = str.find_first_not_of(whitespace);
    if (start == string::npos) {
        return "";
    }
    size_t end = str.find_last_not_of(whitespace);
    return str.substr(start, end - start + 1);
}

double parse_value_with_metric_prefix_util(const string& val_str_orig) {
    string val_str = trim_string_util(val_str_orig);
    if (val_str.empty()) {
        throw invalid_argument("Empty value string for parsing.");
    }
    size_t suffix_start_pos = 0;
    double base_val;
    try {
        base_val = stod(val_str, &suffix_start_pos);
    } catch (const invalid_argument&) {
        throw invalid_argument("Invalid numeric value in string: \"" + val_str_orig + "\"");
    } catch (const out_of_range&) {
        throw out_of_range("Numeric value out of range in string: \"" + val_str_orig + "\"");
    }
    if (suffix_start_pos >= val_str.length()) {
        return base_val;
    }
    string suffix_str = to_lower_util(val_str.substr(suffix_start_pos));
    double multiplier = 1.0;
    if (suffix_str.rfind("meg", 0) == 0) {
        multiplier = 1e6;
    } else if (!suffix_str.empty()) {
        char prefix_char = suffix_str[0];
        switch (prefix_char) {
            case 'p': multiplier = 1e-12; break;
            case 'n': multiplier = 1e-9;  break;
            case 'u': multiplier = 1e-6;  break;
            case 'm': multiplier = 1e-3;  break;
            case 'k': multiplier = 1e3;   break;
            case 'g': multiplier = 1e9;   break;
        }
    }
    return base_val * multiplier;
}

vector<double> gaussian_elimination_matrix(vector<vector<double>> A, vector<double> b) {
    int n = A.size();
    if (n == 0 || (n > 0 && (A[0].size() != n || b.size() != n))) {
        throw runtime_error("Invalid matrix or vector dimensions for Gaussian elimination.");
    }
    for (int i = 0; i < n; ++i) {
        int max_row = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[max_row][i])) max_row = k;
        }
        swap(A[i], A[max_row]);
        swap(b[i], b[max_row]);
        if (abs(A[i][i]) < 1e-12) { }
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[i][i]) < 1e-12) continue;
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) A[k][j] -= factor * A[i][j];
            b[k] -= factor * b[i];
        }
    }
    vector<double> x(n);
    for (int i = n - 1; i >= 0; --i) {
        if (abs(A[i][i]) < 1e-12) {
            if (abs(b[i]) > 1e-12) {
                throw runtime_error("System is inconsistent (no solution).");
            }
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
public:
    string name;
    string node1_name, node2_name;
    double value;
    Component(string name_val, string n1, string n2, double val = 0.0)
            : name(move(name_val)), node1_name(move(n1)), node2_name(move(n2)), value(val) {}
    Component(string name_val, string n1, string n2, const string& val_str)
            : name(move(name_val)), node1_name(move(n1)), node2_name(move(n2)) {
        this->value = parse_value_with_metric_prefix_util(val_str);
    }
    virtual ~Component() = default;
    virtual void stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) = 0;
    virtual string to_netlist_string() const = 0;
};

class VoltageSource;

class Circuit {
public:
    vector<unique_ptr<Component>> components;
    map<string, int> node_to_idx;
    vector<string> idx_to_node_name;
    vector<VoltageSource*> voltage_source_list;
    string ground_node_explicit_name = "0";
    bool ground_node_exists = true;

    void add_component(unique_ptr<Component> comp);
    bool is_ground(const string& node_name) const { return node_name == ground_node_explicit_name; }
    int get_node_matrix_index(const string& node_name) const;
    int get_voltage_source_matrix_index(const string& vs_name);
    int prepare_for_analysis();
    void build_mna_matrix(vector<vector<double>>& A, vector<double>& z);
    void clear();
};

class Resistor : public Component {
public:
    Resistor(const string& r_name, const string& n1, const string& n2, const string& val_str)
            : Component(r_name, n1, n2, val_str) {}
    string to_netlist_string() const override;
    void stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) override;
};

class VoltageSource : public Component {
public:
    VoltageSource(const string& v_name, const string& n1, const string& n2, const string& val_str)
            : Component(v_name, n1, n2, val_str) {}
    string to_netlist_string() const override;
    void stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) override;
};

class CurrentSource : public Component {
public:
    CurrentSource(const string& i_name, const string& n1, const string& n2, const string& val_str)
            : Component(i_name, n1, n2, val_str) {}
    string to_netlist_string() const override;
    void stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) override;
};

void Circuit::clear() {
    components.clear();
    node_to_idx.clear();
    idx_to_node_name.clear();
    voltage_source_list.clear();
    ground_node_exists = true; // Reset to default
    ground_node_explicit_name = "0";
}

void Circuit::add_component(unique_ptr<Component> comp) {
    if (auto vs = dynamic_cast<VoltageSource*>(comp.get())) {
        voltage_source_list.push_back(vs);
    }
    components.push_back(move(comp));
}

int Circuit::get_node_matrix_index(const string& node_name) const {
    if (is_ground(node_name)) return -1;
    auto it = node_to_idx.find(node_name);
    if (it != node_to_idx.end()) {
        return it->second;
    }
    return -2;
}

int Circuit::get_voltage_source_matrix_index(const string& vs_name) {
    for (size_t i = 0; i < voltage_source_list.size(); ++i) {
        if (voltage_source_list[i]->name == vs_name) return i;
    }
    return -1;
}

int Circuit::prepare_for_analysis() {
    node_to_idx.clear();
    idx_to_node_name.clear();
    voltage_source_list.clear();
    set<string> unique_node_names;
    for (const auto& comp : components) {
        unique_node_names.insert(comp->node1_name);
        unique_node_names.insert(comp->node2_name);
        if (auto vs = dynamic_cast<VoltageSource*>(comp.get())) {
            voltage_source_list.push_back(vs);
        }
    }
    int idx = 0;
    for (const auto& name : unique_node_names) {
        if (!is_ground(name)) {
            node_to_idx[name] = idx++;
            idx_to_node_name.push_back(name);
        }
    }
    return idx;
}

void Circuit::build_mna_matrix(vector<vector<double>>& A, vector<double>& z) {
    if (!components.empty() && !ground_node_exists) {
        throw runtime_error("No ground node defined.");
    }
    int N = prepare_for_analysis();
    int M = voltage_source_list.size();
    int system_size = N + M;

    if (system_size == 0) {
        A.clear();
        z.clear();
        return;
    }

    A.assign(system_size, vector<double>(system_size, 0.0));
    z.assign(system_size, 0.0);

    for (const auto& comp : components) {
        comp->stamp(A, z, *this);
    }
}

string Resistor::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}

void Resistor::stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) {
    double conductance = 1.0 / value;
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) A[idx1][idx1] += conductance;
    if (idx2 >= 0) A[idx2][idx2] += conductance;
    if (idx1 >= 0 && idx2 >= 0) {
        A[idx1][idx2] -= conductance;
        A[idx2][idx1] -= conductance;
    }
}

string VoltageSource::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}

void VoltageSource::stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) {
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    int m = circuit.node_to_idx.size() + circuit.get_voltage_source_matrix_index(this->name);
    if (idx1 >= 0) { A[m][idx1] = 1.0; A[idx1][m] = 1.0; }
    if (idx2 >= 0) { A[m][idx2] = -1.0; A[idx2][m] = -1.0; }
    z[m] += this->value;
}

string CurrentSource::to_netlist_string() const {
    ostringstream oss;
    oss << name << " " << node1_name << " " << node2_name << " " << value;
    return oss.str();
}

void CurrentSource::stamp(vector<vector<double>>& A, vector<double>& z, Circuit& circuit) {
    int idx1 = circuit.get_node_matrix_index(node1_name);
    int idx2 = circuit.get_node_matrix_index(node2_name);
    if (idx1 >= 0) z[idx1] -= value;
    if (idx2 >= 0) z[idx2] += value;
}

struct Command {
    string type;
    vector<string> args;
};

class Parser {
public:
    Command parse_line(const string& line);
    void execute_command(Command& cmd, Circuit& circuit);
private:
    void handle_add_component(const vector<string>& args, Circuit& circuit);
    void handle_list_components(const Circuit& circuit);
    void handle_list_nodes(const Circuit& circuit);
};

Command Parser::parse_line(const string& line) {
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

void Parser::execute_command(Command& cmd, Circuit& circuit) {
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
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
}

void Parser::handle_add_component(const vector<string>& args, Circuit& circuit) {
    if (args.size() < 4) throw runtime_error("Insufficient arguments for adding a component. Usage: <name> <node1> <node2> <value>");

    string name = args[0];
    char type_char = toupper(name[0]);
    string n1 = args[1];
    string n2 = args[2];
    string val_str = args[3];

    switch(type_char) {
        case 'R':
            circuit.add_component(make_unique<Resistor>(name, n1, n2, val_str));
            break;
        case 'V':
            circuit.add_component(make_unique<VoltageSource>(name, n1, n2, val_str));
            break;
        case 'I':
            circuit.add_component(make_unique<CurrentSource>(name, n1, n2, val_str));
            break;
        default:
            throw runtime_error("Unknown component type: " + string(1, type_char));
    }
    cout << "Added component " << name << endl;
}

void Parser::handle_list_components(const Circuit& circuit) {
    cout << "Components in circuit:" << endl;
    if (circuit.components.empty()) {
        cout << "  (No components)" << endl;
        return;
    }
    for (const auto& comp : circuit.components) {
        cout << "  - " << comp->to_netlist_string() << endl;
    }
}

void Parser::handle_list_nodes(const Circuit& circuit) {
    set<string> all_nodes;
    for (const auto& comp : circuit.components) {
        all_nodes.insert(comp->node1_name);
        all_nodes.insert(comp->node2_name);
    }
    cout << "Available nodes:" << endl;
    if (all_nodes.empty()) {
        cout << "  (No nodes in circuit)" << endl;
        return;
    }
    for (const auto& node_name : all_nodes) {
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

        line_input = trim_string_util(line_input);
        if (line_input.empty()) continue;

        Command cmd = parser_main.parse_line(line_input);
        parser_main.execute_command(cmd, circuit_main);
    }

    return 0;
}