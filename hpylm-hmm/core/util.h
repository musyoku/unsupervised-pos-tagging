#ifndef _util_
#define _util_
#include <boost/python.hpp>
#include <vector>
#include <iostream>
using namespace std;
using namespace boost;

template<class T>
python::list list_from_vector(vector<T> &vec){  
	 python::list list;
	 typename vector<T>::const_iterator it;

	 for(it = vec.begin(); it != vec.end(); ++it){
		  list.append(*it);
	 }
	 return list;
}

template<class T1,class T2>
python::dict dict_from_map(unordered_map<T1,T2> &map_){  
	 python::dict py_dict;
	 typename unordered_map<T1,T2>::const_iterator it;
	 for(it = map_.begin(); it != map_.end(); ++it){
		  py_dict[it->first]=it->second;        
	 }
	 return py_dict;  
}
double factorial(double n) {
	if (n == 0){
		return 1;
	}
	return n * factorial(n - 1);
}
vector<wstring> split_word_by(const wstring &str, wchar_t delim){
    vector<wstring> elems;
    wstring item;
    for(wchar_t ch: str){
        if (ch == delim){
            if (!item.empty()){
                elems.push_back(item);
            }
            item.clear();
        }
        else{
            item += ch;
        }
    }
    if (!item.empty()){
        elems.push_back(item);
    }
    return elems;
}
void show_progress(int step, int total){
	double progress = step / (double)(total - 1);
	int barWidth = 30;

	cout << "\r" << step << " / " << total << " [";
	int pos = barWidth * progress;
	for(int i = 0; i < barWidth; ++i){
		if (i < pos) cout << "=";
		else if (i == pos) cout << ">";
		else cout << " ";
	}
	cout << "] " << int(progress * 100.0) << " %";
	cout.flush();
	if(step == total){
		cout << endl;
	}
}
#endif