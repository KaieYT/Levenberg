//
//  main.cpp
//  Levenberg
//
//  Created by Yasuda on 2015/06/06.
//  Copyright (c) 2015年 Yasuda. All rights reserved.
//
#include "LM.h"
#include <iostream>

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    LM lm = LM();
    lm.curveFitting("test");
    return 0;
}
