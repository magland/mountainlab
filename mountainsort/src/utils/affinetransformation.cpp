#include "affinetransformation.h"

#include <QList>
#include <math.h>

struct matrix44 {
	double d[4][4];
};

class AffineTransformationPrivate {
public:
	AffineTransformation *q;
	matrix44 m_matrix;

	void set_identity();
	void multiply_by(const matrix44 &M,bool left);
	void copy_from(const AffineTransformation &other);
};

AffineTransformation::AffineTransformation()
{
	d=new AffineTransformationPrivate;
	d->q=this;

	d->set_identity();
}

AffineTransformation::AffineTransformation(const AffineTransformation &other)
{
	d=new AffineTransformationPrivate;
	d->q=this;
	d->copy_from(other);
}

AffineTransformation::~AffineTransformation()
{
	delete d;
}

void AffineTransformation::operator=(const AffineTransformation &other)
{
	d->copy_from(other);
}

Point3D AffineTransformation::map(const Point3D &p)
{
    Point3D q;
	q.x=d->m_matrix.d[0][0]*p.x+d->m_matrix.d[0][1]*p.y+d->m_matrix.d[0][2]*p.z+d->m_matrix.d[0][3];
	q.y=d->m_matrix.d[1][0]*p.x+d->m_matrix.d[1][1]*p.y+d->m_matrix.d[1][2]*p.z+d->m_matrix.d[1][3];
	q.z=d->m_matrix.d[2][0]*p.x+d->m_matrix.d[2][1]*p.y+d->m_matrix.d[2][2]*p.z+d->m_matrix.d[2][3];
	return q;
}

void AffineTransformation::setIdentity()
{
	d->set_identity();
}

void set_identity0(matrix44 &M) {
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			if (i==j) M.d[i][j]=1;
			else M.d[i][j]=0;
		}
	}
}

void copy_matrix0(matrix44 &dst,const matrix44 &src) {
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			dst.d[i][j]=src.d[i][j];
		}
	}
}

void AffineTransformation::rotateX(double theta,bool left)
{
	matrix44 M;
	set_identity0(M);
	M.d[1][1]=cos(theta); M.d[1][2]=-sin(theta);
	M.d[2][1]=sin(theta); M.d[2][2]=cos(theta);
	d->multiply_by(M,left);
}

void AffineTransformation::rotateY(double theta, bool left)
{
	matrix44 M;
	set_identity0(M);
	M.d[0][0]=cos(theta); M.d[0][2]=sin(theta);
	M.d[2][0]=-sin(theta); M.d[2][2]=cos(theta);
	d->multiply_by(M,left);
}

void AffineTransformation::rotateZ(double theta, bool left)
{
	matrix44 M;
	set_identity0(M);
	M.d[0][0]=cos(theta); M.d[0][1]=-sin(theta);
	M.d[1][0]=sin(theta); M.d[1][1]=cos(theta);
	d->multiply_by(M,left);
}

void AffineTransformation::scale(double sx, double sy, double sz, bool left)
{
	matrix44 M;
	set_identity0(M);
	M.d[0][0]=sx;
	M.d[1][1]=sy;
	M.d[2][2]=sz;
	d->multiply_by(M,left);
}

void AffineTransformation::getMatrixData(double *data)
{
	for (int r=0; r<4; r++) {
		for (int c=0; c<4; c++) {
			data[r*4+c]=d->m_matrix.d[r][c];
		}
	}
}

bool AffineTransformation::equals(const AffineTransformation &other)
{
	for (int r=0; r<4; r++) {
		for (int c=0; c<4; c++) {
			if (d->m_matrix.d[r][c]!=other.d->m_matrix.d[r][c]) return false;
		}
	}
	return true;
}



QList<double> invert33(QList<double> &data33) {
	double X1[3][3];
	double X2[3][3];
	int ct=0;
	for (int j=0; j<3; j++)
	for (int i=0; i<3; i++) {
		X1[i][j]=data33[ct];
		ct++;
	}

	X2[0][0] = X1[1][1]*X1[2][2] - X1[2][1]*X1[1][2];
	X2[0][1] = X1[0][1]*X1[2][2] - X1[2][1]*X1[0][2];
	X2[0][2] = X1[0][1]*X1[1][2] - X1[1][1]*X1[0][2];
	X2[1][0] = X1[1][0]*X1[2][2] - X1[2][0]*X1[1][2];
	X2[1][1] = X1[0][0]*X1[2][2] - X1[2][0]*X1[0][2];
	X2[1][2] = X1[0][0]*X1[1][2] - X1[1][0]*X1[0][2];
	X2[2][0] = X1[1][0]*X1[2][1] - X1[2][0]*X1[1][1];
	X2[2][1] = X1[0][0]*X1[2][1] - X1[2][0]*X1[0][1];
	X2[2][2] = X1[0][0]*X1[1][1] - X1[1][0]*X1[0][1];

	X2[1][0]*=-1; X2[0][1]*=-1; X2[2][1]*=-1; X2[1][2]*=-1;

	double det= 	 X1[0][0]*X1[1][1]*X1[2][2]
				+X1[0][1]*X1[1][2]*X1[2][0]
				+X1[0][2]*X1[1][0]*X1[2][1]
				-X1[0][2]*X1[1][1]*X1[2][0]
				-X1[0][0]*X1[1][2]*X1[2][1]
				-X1[0][1]*X1[1][0]*X1[2][2] ;
	if (det!=0) {
		for (int j=0; j<3; j++)
		for (int i=0; i<3; i++)
			X2[i][j]/=det;
	}

	QList<double> ret;
	for (int j=0; j<3; j++)
	for (int i=0; i<3; i++) {
		ret << X2[i][j];
	}
	return ret;
}

AffineTransformation AffineTransformation::inverse() const
{
	AffineTransformation T;

	QList<double> data33;
	for (int j=0; j<3; j++)
	for (int i=0; i<3; i++)
		data33 << d->m_matrix.d[i][j];
	data33=invert33(data33);

	double tmp[4][4];
	int ct=0;
	for (int j=0; j<3; j++)
	for (int i=0; i<3; i++) {
		tmp[i][j]=data33[ct];
		ct++;
	}
	tmp[0][3]=-(tmp[0][0]*d->m_matrix.d[0][3]+tmp[0][1]*d->m_matrix.d[1][3]+tmp[0][2]*d->m_matrix.d[2][3]);
	tmp[1][3]=-(tmp[1][0]*d->m_matrix.d[0][3]+tmp[1][1]*d->m_matrix.d[1][3]+tmp[1][2]*d->m_matrix.d[2][3]);
	tmp[2][3]=-(tmp[2][0]*d->m_matrix.d[0][3]+tmp[2][1]*d->m_matrix.d[1][3]+tmp[2][2]*d->m_matrix.d[2][3]);
	tmp[3][0]=0; tmp[3][1]=0; tmp[3][2]=0; tmp[3][3]=1;
	QList<double> data44;
	for (int i=0; i<4; i++)
	for (int j=0; j<4; j++) {
		T.d->m_matrix.d[i][j]=tmp[i][j];
	}

	return T;
}

void AffineTransformationPrivate::set_identity()
{
	set_identity0(m_matrix);
}

void AffineTransformationPrivate::multiply_by(const matrix44 &M, bool left)
{
	matrix44 tmp; copy_matrix0(tmp,m_matrix);
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			double val=0;
			if (left) {
				for (int k=0; k<4; k++) val+=M.d[i][k]*tmp.d[k][j];
			}
			else {
				for (int k=0; k<4; k++) val+=tmp.d[i][k]*M.d[k][j];
			}
			m_matrix.d[i][j]=val;
		}
	}
}

void AffineTransformationPrivate::copy_from(const AffineTransformation &other)
{
	for (int j=0; j<4; j++) {
		for (int i=0; i<4; i++) {
			m_matrix.d[i][j]=other.d->m_matrix.d[i][j];
		}
	}
}
