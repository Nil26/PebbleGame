//
// Created by Shang Zhang on 4/13/17.
//

#ifndef PEBBLEGAMETEST_BOND_H
#define PEBBLEGAMETEST_BOND_H

#include <iostream>

class Bond{
public:
    std::pair<int, int> vertices; // the connected vertices
    int RigidIndex;         //  the rigid cluster index
    Bond (int x_input, int y_input);
    void initBondRigidIndex ();
    std::pair<int,int> vertex() const { return vertices;}; // to get the vertices of the class Bond
};

#endif //PEBBLEGAMETEST_BOND_H
