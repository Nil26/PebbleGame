//
//  main.cpp
//  SiteRP_cont
//
//  Created by Shang Zhang on 5/8/17.
//  Copyright © 2017 Shang Zhang. All rights reserved.
//

#include <iostream>
#include <cstring>
#include <cstdlib>
#include "SiteRP_cont.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    if (argc < 4) { // We expect 5 arguments: the program name, number of particles, source file name, the source path and the destination path
        std::cerr << "Usage: " << argv[0] << " <PARTICLE NUMBER> <SOURCE NAME> <SOURCE PATH> <DESTINATION PATH> <BOND CRITERIA-(OPTIONAL)>" << std::endl;
        return 1;
    }
    
    int NumParticle = atoi(argv[1]);
    std::string FileName = std::string(argv[2]);
    std::string ReadPATH = std::string(argv[3]);
    std::string OutPATH = std::string(argv[4]);
    double BondLEN = atof(argv[5]);
//    std::string ReadPATH = "./MCdata/DataQuenchT0.2/B0QuenchTin0_5to0_2dt1e6/N10000/";
//    std::string OutPATH = "./data/DataQuenchT0.2/B0QuenchTin0_5to0_2dt1e6/N10000/";
//    std::string FileName = "cfgdensity0_48";
//    int NumParticle = 10000;
    //SiteRP_cont(int sizeIN, std::string FileNameIN, std::string ReadPATHIN, std::string OutPATHIN)
    SiteRP_cont a(NumParticle, FileName, ReadPATH, OutPATH, BondLEN);
    
    //a.ReadVerticeInfoFromCSV();
    //a.StudyArtificialNetwork();
    a.ReadVerticeInfofromdump();
    //a.ReadVerticeInfo();
    a.BuildNetwork();
    //a.CoordNumber();
}
