/******************************************************
** See the accompanying README and LICENSE files
** Author(s): Jeremy Magland
*******************************************************/

#ifndef AFFINETRANSFORMATION_H
#define AFFINETRANSFORMATION_H

struct Point3D {
    Point3D(float x0 = 0, float y0 = 0, float z0 = 0)
    {
        x = x0;
        y = y0;
        z = z0;
    }

    float x, y, z;
};

class AffineTransformationPrivate;
class AffineTransformation {
public:
    friend class AffineTransformationPrivate;
    AffineTransformation();
    AffineTransformation(const AffineTransformation& other);
    ~AffineTransformation();
    void operator=(const AffineTransformation& other);

    Point3D map(const Point3D& p);

    void setIdentity();
    void rotateX(double theta, bool left = true);
    void rotateY(double theta, bool left = true);
    void rotateZ(double theta, bool left = true);

    void scale(double sx, double sy, double sz, bool left = true);
    void getMatrixData(double* data);

    bool equals(const AffineTransformation& other);

    AffineTransformation inverse() const;

private:
    AffineTransformationPrivate* d;
};

#endif // AFFINETRANSFORMATION_H
