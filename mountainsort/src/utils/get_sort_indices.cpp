#include "get_sort_indices.h"

#include <QVector>

QList<long> get_sort_indices(const QList<long>& X)
{
    QList<long> result;
    result.reserve(X.size());
    for (int i = 0; i < X.size(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
                     [&X](int i1, int i2) { return X[i1] < X[i2]; });
    return result;
}

QList<long> get_sort_indices(const QVector<double>& X)
{
    QList<long> result;
    result.reserve(X.size());
    for (int i = 0; i < X.size(); ++i)
        result << i;
    std::stable_sort(result.begin(), result.end(),
                     [&X](int i1, int i2) { return X[i1] < X[i2]; });
    return result;
}
