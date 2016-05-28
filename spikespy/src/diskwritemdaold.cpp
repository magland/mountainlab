#include "diskwritemdaold.h"
#include <QDir>
#include <QDebug>
#include "usagetracking.h"
#include "sscommon.h"
#include "sscontroller.h"

class DiskWriteMdaOldPrivate {
public:
    DiskWriteMdaOld *q;
	QString m_path;
	int m_data_type;
	int m_total_size;
	bool m_allocated;
	FILE *m_file;
	MDAIO_HEADER m_mda_header;

	void clear_all();
	int get_index(int i1,int i2,int i3=0,int i4=0,int i5=0,int i6=0);
};


DiskWriteMdaOld::DiskWriteMdaOld(const QString &path)
{
    d=new DiskWriteMdaOldPrivate;
	d->q=this;
	d->m_data_type=MDAIO_TYPE_FLOAT32;
	d->m_allocated=false;
	d->m_file=0;
	d->m_total_size=1;

	if (!path.isEmpty()) setPath(path);
	else useTemporaryFile();
}

DiskWriteMdaOld::~DiskWriteMdaOld()
{
	d->clear_all();
	if (d->m_file) jfclose(d->m_file);
	delete d;
}

void DiskWriteMdaOld::setPath(const QString &path)
{
	if (d->m_file) jfclose(d->m_file);
	d->m_file=0;
	d->m_path=path;
	d->clear_all();
}

void DiskWriteMdaOld::useTemporaryFile()
{
    QString path=QString(ssTempPath()+"/spikespy.diskwritemdaold.tmp.%1").arg(qAbs(qrand()));
	removeOnClose(path);
	setPath(path);
}

void DiskWriteMdaOld::setDataType(int data_type)
{
	if (d->m_allocated) {
		printf("Warning: cannot set data type once array has been allocated.\n");
		return;
	}
	d->m_data_type=data_type;
}

void DiskWriteMdaOld::allocate(int N1, int N2, int N3, int N4, int N5, int N6)
{
	if (d->m_allocated) {
		printf("Warning: cannot allocate once array has been allocated.\n");
		return;
	}
	jfclose(jfopen(d->m_path.toLatin1().data(),"wb")); //remove the file
	d->m_file=jfopen(d->m_path.toLatin1().data(),"wb+");
	if (!d->m_file) {
		printf("Warning: cannot allocate. Unable to open file: %s\n",d->m_path.toLatin1().data());
		return;
	}
	int dims[MDAIO_MAX_DIMS]; for (int i=0; i<MDAIO_MAX_DIMS; i++) dims[i]=1;
	dims[0]=N1;
	dims[1]=N2;
	dims[2]=N3;
	dims[3]=N4;
	dims[4]=N5;
	dims[5]=N6;
	int num_dims=2;
	for (int i=2; i<MDAIO_MAX_DIMS; i++) {
		if (dims[i]>1) num_dims=i+1;
	}
	d->m_total_size=1;
	for (int i=0; i<num_dims; i++) {
		d->m_total_size*=dims[i];
	}
	d->m_mda_header.data_type=d->m_data_type;
	for (int i=0; i<MDAIO_MAX_DIMS; i++) d->m_mda_header.dims[i]=dims[i];
	d->m_mda_header.num_dims=num_dims;
	mda_write_header(&d->m_mda_header,d->m_file);
    double *all_zeros=(double *)jmalloc(sizeof(double)*d->m_total_size);
	for (int j=0; j<d->m_total_size; j++) all_zeros[j]=0;
    mda_write_float64(all_zeros,&d->m_mda_header,d->m_total_size,d->m_file);
	jfree(all_zeros);
	d->m_allocated=true;
}

void DiskWriteMdaOld::setValue(double val, int i1, int i2, int i3, int i4, int i5,int i6)
{
	if (!d->m_file) return;
	int ind=d->get_index(i1,i2,i3,i4,i5,i6);
	if ((ind>=0)&&(ind<d->m_total_size)) {
		fseek(d->m_file,d->m_mda_header.header_size+ind*d->m_mda_header.num_bytes_per_entry,SEEK_SET);
		if (d->m_data_type==MDAIO_TYPE_BYTE) {
			unsigned char tmp=(unsigned char)val;
			fwrite(&tmp,d->m_mda_header.num_bytes_per_entry,1,d->m_file);
		}
		else if (d->m_data_type==MDAIO_TYPE_FLOAT32) {
			float tmp=(float)val;
			fwrite(&tmp,d->m_mda_header.num_bytes_per_entry,1,d->m_file);
		}
		else if (d->m_data_type==MDAIO_TYPE_INT16) {
			int16_t tmp=(int16_t)val;
			fwrite(&tmp,d->m_mda_header.num_bytes_per_entry,1,d->m_file);
		}
		else if (d->m_data_type==MDAIO_TYPE_INT32) {
			int32_t tmp=(int32_t)val;
			fwrite(&tmp,d->m_mda_header.num_bytes_per_entry,1,d->m_file);
		}
		else if (d->m_data_type==MDAIO_TYPE_UINT16) {
			quint16 tmp=(quint16)val;
			fwrite(&tmp,d->m_mda_header.num_bytes_per_entry,1,d->m_file);
		}
        else if (d->m_data_type==MDAIO_TYPE_FLOAT64) {
            double tmp=(double)val;
            fwrite(&tmp,d->m_mda_header.num_bytes_per_entry,1,d->m_file);
        }
	}
}

DiskReadMdaOld DiskWriteMdaOld::toReadMda()
{
    return DiskReadMdaOld(d->m_path);
}



void DiskWriteMdaOldPrivate::clear_all()
{
	if (m_file) jfclose(m_file);
	m_file=0;
}

int DiskWriteMdaOldPrivate::get_index(int i1, int i2, int i3, int i4, int i5, int i6)
{
	int inds[6]; inds[0]=i1; inds[1]=i2; inds[2]=i3; inds[3]=i4; inds[4]=i5; inds[5]=i6;
	int factor=1;
	int ret=0;
	for (int j=0; j<6; j++) {
		ret+=factor*inds[j];
		factor*=m_mda_header.dims[j];
	}
	return ret;
}
