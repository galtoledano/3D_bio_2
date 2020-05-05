//
// Created by galtoledano on 05/05/2020.
//

#include "structalign.h"
#include "Ex2Code/Vector3.h"
#include "Ex2Code/Atom.h"
#include "Ex2Code/RigidTrans3.h"
#include "Ex2Code/Matrix3.h"
#include "Ex2Code/Molecule.h"
#include "Ex2Code/PDB.h"
#include "Ex2Code/Match.h"
#include "Ex2Code/GeomHash.h"
#include <chrono>
#include <iostream>

void RigidTrans3(const Vector3 &a, const Vector3 &b, const Vector3 &c)
{
    Vector3 translation = (a+b+c)/3.0;
    Vector3 x = translation - a;
    x /= x.norm();
    Vector3 z = (c-a)&(c-b);
    z /= z.norm();
    Vector3 y = z & x;
    Matrix3 rotation = Matrix3(x,y,z);
    rotation.transpose();
}

int main(int argc , char* argv[]){
    // measure the run time
    auto start = std::chrono::system_clock::now();

    if(argc != 4) {
        std::cerr << "Usage: "<<argv[0]<< " target.pdb model.pdb num_rotations dist_threshold" << std::endl;
        exit(1);
    }

    //********Parameters********************
    float m_fDistThr = atof(argv[1]); // distance threshold on atoms in correspondence

    std::cout << "Distance threshold: "<< m_fDistThr  << std::endl;

    // read the two files into Molecule
    Molecule<Atom> molModel, molTarget;

    std::ifstream fileModel(argv[2]);
    std::ifstream fileTarget(argv[3]);

    if(!fileModel) {
        std::cout<< "File " << argv[2] << " does not exist." << std::endl;
        return 0;
    }
    if(!fileTarget) {
        std::cout << "File " << argv[3] << " does not exist." << std::endl;
        return 0;
    }

//    molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
//    molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());

//    // calculate center of mass
//    Vector3 vectModelMass(0,0,0);
//    for(unsigned int i=0; i <molModel.size(); i++) {
//        vectModelMass+=molModel[i].position();
//    }
//    vectModelMass/=molModel.size();
//
//    Vector3 vectTargetMass(0,0,0);
//    for(unsigned int i=0; i<molTarget.size(); i++) {
//        vectTargetMass+=molTarget[i].position();
//    }
//    vectTargetMass/=molTarget.size();
//
//    // transform the molecules to the center of the coordinate system
//    molModel+=(-vectModelMass);
//    molTarget+=(-vectTargetMass);
//
//    // next we insert the target molecule into hash
//    // this will help us to find atoms that are close faster
//    GeomHash <Vector3,int> gHash(3,m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
//    for(unsigned int i=0; i<molTarget.size(); i++) {
//        gHash.insert(molTarget[i].position(), i); // coordinate is the key to the hash, we store atom index
//    }
//
    printf("hello ! \n");


}
