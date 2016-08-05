#include "get_sort_indices.h"

#include <QVector>
#include <algorithm>

QList<long> get_sort_indices(const QList<long>& X)
{
    QList<long> result = X;
    std::stable_sort(result.begin(), result.end());
    return result;
}

QList<long> get_sort_indices(const QVector<double>& X)
{
    QList<long> result;
    std::copy(X.constBegin(), X.constEnd(), result.begin());
    std::stable_sort(result.begin(), result.end());
    return result;
}
