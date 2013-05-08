//
//  Parser.h
//  PHOTONVIZ
//
//  Created by Jaromir Herskovic on 06.05.13.
//  Copyright (c) 2013 ZDENEK GLAZER. All rights reserved.
//

#ifndef __PHOTONVIZ__Parser__
#define __PHOTONVIZ__Parser__

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <istream>
#include <sstream>
#include <iostream>
#include <fstream> 
#include <sstream>
#include <string>
#include "mystructs.h"

class Parser {
    template<typename T>
    T StringToNumber(const std::string& numberAsString)
    {
        T valor;
        
        std::stringstream stream(numberAsString);
        stream >> valor;
        if (stream.fail()) {
            std::runtime_error e(numberAsString);
            throw e;
        }
        return valor;
    }

public:
    Parser(){};
    bool getLines(std::vector<volumePath*> &lns,  std::string fn){
        std::ifstream file(fn.c_str() );
        
        std::string line;
        while( std::getline( file, line ) )
        {
            std::istringstream iss( line );
            
            std::string result;
            if( std::getline( iss, result , '=') )
            {
                if( result == "l" )
                {
                    volumePath* vpth= new volumePath();
                    std::string token;
                    bool cond=true;
                    float* point;
                    while (cond) {
                        point=new float[4];
                        for (int i=0; i<4; ++i) {
                            if(!std::getline( iss, token, ',' )){
                                delete[] point;
                                point=NULL;
                                cond=false;
                                break;
                            }
                            point[i]=StringToNumber<float>(token);
                        }
                        if (point==NULL) {
                            break;
                        }
                        vpth->points.push_back(point);
                    }
                    lns.push_back(vpth);
                }
            }
        }
        return true;
    }
    
    bool getPoints(std::vector<float*> &points, std::string fileName){
        std::ifstream file(fileName.c_str() );
        
        std::string line;
        while( std::getline( file, line ) )
        {
            std::istringstream iss( line );
            
            std::string result;
            if( std::getline( iss, result , '=') )
            {
                if( result == "p" )
                {
                    float* point=new float[3];
                    std::string token;
                    for (int i=0; i<3; ++i) {
                        if(!std::getline( iss, token, ',' )){
                            delete[] point;
                            break;
                        }
                        point[i]=StringToNumber<float>(token);
                    }
                    if (point==NULL) {
                        continue;
                    }
                    points.push_back(point);
                }
            }
        }
        return true;
    }
};



#endif /* defined(__PHOTONVIZ__Parser__) */
