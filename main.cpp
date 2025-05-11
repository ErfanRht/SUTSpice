#include <iostream>
#include <string>
#include <vector>
#include <memory>
using namespace std;

class Circuit;

class Component {
public:
    string name;
    string node1_name, node2_name;
    double value;

    Component(string name_val, string n1, string n2, double val = 0.0)
            : name(move(name_val)), node1_name(move(n1)), node2_name(move(n2)), value(val) {}

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
    Circuit circuit;
    cout << "SUTSpice - Basic classes created." << endl;
    cout << "A circuit has been created." << endl;
    return 0;
}