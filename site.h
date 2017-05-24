//
// Created by Shang Zhang on 4/20/17.
//

#ifndef PEBBLEGAMETEST_SITE_H
#define PEBBLEGAMETEST_SITE_H

class site{
public:
    double XCoordinate;
    double YCoordinate;
    int SiteIndex;  // the index from Monte-Carlo data, not in the vector vertices!!!

    site(double xCoord, double yCoord) {XCoordinate = xCoord; YCoordinate = yCoord;};
    site(double xCoord, double yCoord, int Index) {XCoordinate = xCoord; YCoordinate = yCoord; SiteIndex = Index;};
};

#endif //PEBBLEGAMETEST_SITE_H
