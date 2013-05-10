/*
 *  datadumper.h
 *  pbrt
 *
 *  Created by Zdenek Glazer on 12/1/12.
 *  Copyright 2012 CVUT. All rights reserved.
 *
 */
#ifndef DATA_DUMPER_H
#define DATA_DUMPER_H
#include "igi.h"
#include "pbrt.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <ostream>
//#include "vlstructs.h"
struct VirtualSphericalLight;
struct VolumePath;
inline bool toLineString(VolumePath &vpth,string &data);
inline string toPointString(VirtualSphericalLight &vsl);

using namespace std;
class DataDumper {
	public:
	DataDumper(){}
	~DataDumper(){}
    void dumpVirtualSphericalLight(std::vector<VirtualSphericalLight*> &vsls,const string &fileName){
        ofstream myfile;
		printf("data size is %d",vsls.size());
		myfile.open(fileName.c_str());
		VirtualSphericalLight* curl;
		for (int i=0; i<vsls.size(); i++) {
            curl=vsls[i];
            myfile<<"p=";
			myfile<<toPointString(*curl);
			myfile<<"\n";
		}
		//myfile << "Writing this to a file.\n";
		myfile.close();
    }
    
    void dumpVolumePaths(std::vector<VolumePath*>vpths,const string &fileName){
        ofstream myfile;
		printf("data size is %d",vpths.size());
		myfile.open(fileName.c_str());
        VolumePath* curl;
        std::string str;
		for (int i=0; i<vpths.size(); i++) {
            curl=vpths[i];
            if(!toLineString(*curl,str)) continue;
            myfile<<"l=";
			myfile<<str;
			myfile<<"\n";
		}
		//myfile << "Writing this to a file.\n";
		myfile.close();

    }
	void dump2File(const string &fileName){
		ofstream myfile;
		printf("data size is %d",data.size());
		myfile.open(fileName.c_str());
		myfile<<"{\n\"data\":{\n";
		for (it=data.begin(); it!=data.end(); it++) {
			myfile<<"\"";
			myfile<<it->first;
			myfile<<"\"";
			myfile<<":[";
			myfile<<it->second;
			myfile<<"]";
			myfile<<"\n";
		}
		myfile<<"\n}\n}";
		//myfile << "Writing this to a file.\n";
		myfile.close();
	}
	void dump(string keyString,string dataString){
		// there is no such data so i will have to start one
		it=data.find(keyString);
		if (it==data.end()) {
			data.insert(pair<string, string>(keyString, dataString));
			//printf("new data type: %s\n",keyString.c_str());
		}else {
			it->second=it->second+","+dataString;
			//printf("inserting [%s,%s] ",keyString.c_str(),dataString.c_str());
		}

	}
	string vpl2String();
private:
	int shit;
	map <string, string> data;
	map<string,string>::iterator it;
};

#endif