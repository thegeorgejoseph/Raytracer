#ifndef VECTOR_H
#define VECTOR_H

struct Vector {
    double x, y, z;        

    Vector(double a = 0, double b = 0, double c = 0){
        x = a;
        y = b;
        z = c;
    }

    Vector operator+(const Vector &b) const {
        return Vector(x + b.x, y + b.y, z + b.z);
    }

    Vector operator*(double b) const { 
        return Vector(x * b,y * b,z * b); 
    }

    Vector operator/(double b) const { 
        return Vector(x / b,y / b,z / b); 
    }

    Vector operator-(const Vector &b) const { 
        return Vector(x - b.x, y - b.y, z - b.z); 
    }

    Vector operator-() const { 
        return Vector(-x,-y,-z); 
    }
};

#endif