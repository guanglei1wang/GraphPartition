//
// Created by Hassan on 19/11/2015.
//

#ifndef GRAVITY_CONSTANT_H
#define GRAVITY_CONSTANT_H
#include <iostream>
#include <vector>
#include <forward_list>
#include <assert.h>
#include <string>
#include <map>
#include <memory>
#include <gravity/types.h>


using namespace std;



/** Backbone class for constant */
class constant_{
protected:
    CType               _type;
    
    
public:
    bool                _is_transposed = false; /**< True if the constant is considered as a transposed vector */
    bool                _is_vector = false; /**< True if the constant is considered as a vector */
    size_t              _dim = 0; /*<< dimension of current vector */
    
    virtual ~constant_(){};
    CType get_type() const { return _type;}
    void set_type(CType type){ _type = type;}
    
    /** Querries */
    bool is_binary() const{
        return (_type==binary_c);
    };
    
    bool is_short() const{
        return (_type==short_c);
    }
    
    bool is_integer() const{
        return (_type==integer_c);
    };
    
    bool is_float() const{
        return (_type==float_c);
    };
    
    bool is_double() const{
        return (_type==double_c);
    };
    
    bool is_long() const{
        return (_type==long_c);
    };
    
    bool is_number() const{
        return (_type!=par_c && _type!=uexp_c && _type!=bexp_c && _type!=var_c && _type!=sdpvar_c && _type!=func_c);
    }
    bool is_param() const{
        return (_type==par_c);
    };

    bool is_uexpr() const{
        return (_type==uexp_c);
    };

    bool is_bexpr() const{
        return (_type==bexp_c);
    };
    
    bool is_expr() const{
        return (_type==uexp_c || _type==bexp_c);
    };

    
    bool is_var() const{
        return (_type==var_c);
    };
    
    bool is_sdpvar() const{
        return (_type==sdpvar_c);
    };
    
    bool is_function() const{
        return (_type==func_c);
    };
    
    
    Sign get_all_sign() const;
    Sign get_sign(int idx=0) const;
    bool is_zero() const; /**< Returns true if constant equals 0 */
    bool is_unit() const; /**< Returns true if constant equals 1 */
    bool is_neg_unit() const; /**< Returns true if constant equals -1 */
    bool is_positive() const; /**< Returns true if constant is positive */
    bool is_negative() const; /**< Returns true if constant is negative */
    bool is_non_positive() const; /**< Returns true if constant is non positive */
    bool is_non_negative() const; /**< Returns true if constant is non negative */
};

template<typename type>
class param;

/** Polymorphic class constant, can store an arithmetic number (int. float, double..).*/
template<typename type = float>
class constant: public constant_{
protected:
    type        _val;
public:
    
    /** Constructors */
    constant(){        
        if(typeid(type)==typeid(bool)){
            set_type(binary_c);
            return;
        }
        if(typeid(type)==typeid(short)) {
            set_type(short_c);
            return;
        }
        if(typeid(type)==typeid(int)) {
            set_type(integer_c);
            return;
        }
        if(typeid(type)==typeid(float)) {
            set_type(float_c);
            return;
        }
        if(typeid(type)==typeid(double)) {
            set_type(double_c);
            return;
        }
        if(typeid(type)==typeid(long double)) {
            set_type(long_c);
            return;
        }
        throw invalid_argument("Unknown constant type.");
    }
    
    constant(const constant& c){ /**< Copy constructor */
        _type = c._type;
        _val = c._val;
        _is_transposed = c._is_transposed;
        _is_vector = c._is_vector;
    };

    constant(type val):constant(){
        _val = val;
    };
    


    ~constant(){};
    
    constant tr(){
        auto newc(*this);
        newc._is_transposed = true;
        return newc;
    };
    
    type eval() const { return _val;}
    
    void set_val(type val) {
        _val = val;
    }
    
    Sign get_sign() const{
        if (_val==0) {
            return zero_;
        }
        if (_val > 0) {
            return pos_;
        }
        if (_val < 0) {
            return neg_;
        }
        return unknown_;
    }
    
    /** Operators */
    bool is_negative() const {
        return _val < 0;
    }
    
    bool operator==(const constant& c) const {
        return (_type==c._type && _val==c._val);
    }
    
    bool operator==(const type& v) const{
        return _val==v;
    }

    constant& operator=(const type& val){
        _val = val;
        return *this;
    }
    
    constant& operator+=(const type& v){
        _val += v;
        return *this;
    }
    
    constant& operator-=(const type& v){
        _val -= v;
        return *this;
    }
    
    constant& operator*=(const type& v){
        _val *= v;
        return *this;
    }
    
    constant& operator/=(const type& v){
        _val /= v;
        return *this;
    }
    
    friend constant operator+(const constant& c1, const constant& c2){
        return constant(c1._val + c2._val);
    }

    friend constant operator-(const constant& c1, const constant& c2){
        return constant(c1._val - c2._val);
    }
    
    friend constant operator/(const constant& c1, const constant& c2){
        return constant(c1._val / c2._val);
    }
    
    friend constant operator*(const constant& c1, const constant& c2){
        return constant(c1._val * c2._val);
    }
    
    friend constant operator^(const constant& c1, const constant& c2){
        return constant(pow(c1._val,c2._val));
    }


    friend constant operator+(const constant& c, type cst){
        return constant(c._val + cst);
    }
    
    friend constant operator-(const constant& c, type cst){
        return constant(c._val - cst);
    }
    
    friend constant operator*(const constant& c, type cst){
        return constant(c._val * cst);
    }

    
    friend constant operator/(const constant& c, type cst){
        return constant(c._val / cst);
    }

    friend constant operator+(type cst, const constant& c){
        return constant(c._val + cst);
    }
    
    friend constant operator-(type cst, const constant& c){
        return constant(cst - c._val);
    }
    
    friend constant operator*(type cst, const constant& c){
        return constant(c._val * cst);
    }
    
    
    friend constant operator/(type cst, const constant& c){
        return constant(cst / c._val);
    }

    friend constant cos(const constant& c){
        return constant(cos(c._val));
    }
    
    friend constant sin(const constant& c){
        return constant(sin(c._val));
    }
    
    friend constant sqrt(const constant& c){
        return constant(sqrt(c._val));
    }
    
    friend constant expo(const constant& c){
        return constant(exp(c._val));
    }
    
    friend constant log(const constant& c){
        return constant(log(c._val));
    }

    
    /** Output */
    void print() const{
        cout << _val;
    }
    
    string to_str() const{
        char buffer [50];
        if(typeid(type)==typeid(float) || typeid(type)==typeid(double) || typeid(type)==typeid(long double)){
            sprintf (buffer, "%g", _val);
        }
        else {
         sprintf (buffer, "%d", _val);
        }
//        cout << string(buffer) << endl;
        return string(buffer);
//        return std::to_string(_val);
    }

    
};



#endif //GRAVITY_CONSTANT_H
