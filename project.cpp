#pragma once 
#include <TBRRT.h>


int main(){

    //set the start position in joint space

    //set the end position in joint space
    //create an instance of the rrt object
    int maxNodes = 10000;

    TBRRT rrt(start,end);

    path = rrt.search(maxNodes);


}