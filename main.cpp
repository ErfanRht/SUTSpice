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
    virtual void stamp(vector<vector<double>>& G, Circuit& circuit) = 0;
    virtual string to_netlist_string() const = 0;
};

class Circuit {
public:
    vector<unique_ptr<Component>> components;
    map<string, int> node_to_idx;
    vector<string> idx_to_node_name;
    bool ground_node_exists = false;
    string ground_node_explicit_name = "0";

    void add_component(unique_ptr<Component> comp) {
        components.push_back(move(comp));
    }

    bool is_ground(const string& node_name) const {
        return node_name == ground_node_explicit_name;
    }

    int get_node_matrix_index(const string& node_name) {
        if (is_ground(node_name)) return -1;
        if (node_to_idx.find(node_name) == node_to_idx.end()) {
            int new_idx = node_to_idx.size();
            node_to_idx[node_name] = new_idx;
            idx_to_node_name.push_back(node_name);
        }
        return node_to_idx.at(node_name);
    }
};

class Resistor : public Component {
public:
    Resistor(const string& r_name, const string& n1, const string& n2, const string& val_str)
            : Component(r_name, n1, n2, val_str) {
        if (value <= 0) throw runtime_error("Resistor " + name + " must have positive resistance.");
    }

    string to_netlist_string() const override {
        ostringstream oss;
        oss << name << " " << node1_name << " " << node2_name << " " << fixed << setprecision(12) << value;
        return oss.str();
    }

    void stamp(vector<vector<double>>& G, Circuit& circuit) override {
        double conductance = 1.0 / value;
        int idx1 = circuit.get_node_matrix_index(node1_name);
        int idx2 = circuit.get_node_matrix_index(node2_name);

        if (idx1 >= 0) G[idx1][idx1] += conductance;
        if (idx2 >= 0) G[idx2][idx2] += conductance;
        if (idx1 >= 0 && idx2 >= 0) {
            G[idx1][idx2] -= conductance;
            G[idx2][idx1] -= conductance;
        }
    }
};

int main() {
    Circuit circuit;
    circuit.add_component(make_unique<Resistor>("R1", "1", "0", "1k"));
    cout << "SUTSpice - Resistor implemented with stamping." << endl;
    cout << circuit.components[0]->to_netlist_string() << endl;
    return 0;
}