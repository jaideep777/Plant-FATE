#ifndef TIME_ARITHMETICS_H
#define TIME_ARITHMETICS_H

#include <iostream>
#include <boost/date_time/gregorian/gregorian.hpp>

using namespace boost::gregorian;


double year_diff(date d1, date d2){

    date_period dp(d1, d2);
    std::cout << dp << std::endl;
    double a = (double) dp.length().days();

    int y1 = d1.year();
    std::cout << y1 << std::endl;
    bool y1_isleapyr = gregorian_calendar::is_leap_year(y1);
    std::cout << y1_isleapyr << std::endl;
    int y2 = d1.year();
    std::cout << y2 << std::endl;
    bool y2_isleapyr = gregorian_calendar::is_leap_year(y2);
    std::cout << y2_isleapyr << std::endl;
    
    double yr_diff = 0;
    if(y1 == y2){
        yr_diff = a / 365;
        std::cout << yr_diff<< std::endl;
        if(y1_isleapyr){
            yr_diff = a /366;
        }
        return yr_diff;
    }    

    double yr1_diff = date_period(d1, date(y1, 12, 31)).length().days();
    if(y1_isleapyr){
        yr1_diff /= 366;
    }
    else{
        yr1_diff /=365;
    }
    double yr2_diff = date_period(date(y2, 1, 1), d2).length().days();
    if(y2_isleapyr){
        yr2_diff /= 366;
    }
    else{
        yr2_diff /=365;
    }

    yr_diff = yr1_diff + yr2_diff + (y2 - y1 - 1);
    return yr_diff;
}


#endif