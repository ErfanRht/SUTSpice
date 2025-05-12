#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <algorithm>
#include <cctype>
#include <stdexcept>
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
};

class Circuit {
public:
    vector<unique_ptr<Component>> components;

    void add_component(unique_ptr<Component> comp) {
        components.push_back(move(comp));
    }
};

int main() {
    cout << "SUTSpice - Utility functions added." << endl;
    try {
        double val = parse_value_with_metric_prefix_util("1.5k");
        cout << "Parsed '1.5k' to: " << val << endl;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
    }
    return 0;
}