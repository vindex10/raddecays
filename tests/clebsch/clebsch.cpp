#include <iostream>
#include "utils.hpp"

using namespace std;

int main() {
    cout << clebsch(2.*1.+1., 2.*0.+1., 2.*1.+1., 2.*1.+1., 2.*2.+1., 2.*1.+1.) << endl;
    return 0;
}
