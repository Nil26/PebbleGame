//
// Created by Shang Zhang on 4/13/17.
//

#include "bond.h"

void Bond::initBondRigidIndex ()
{
    RigidIndex = 0;
}

Bond::Bond (int x_input, int y_input)
{
    vertices = std::make_pair(x_input,y_input);
    RigidIndex = 0;
}