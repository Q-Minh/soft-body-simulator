#include <sbs/math/sph.h>

int main(int argc, char** argv)
{
    Eigen::Vector3d Xi, X1, X2, X3;

    Xi << 0., 0., 0.;
    X1 << 1., 0., 0.;
    X2 << 0., 1., 0.;
    X3 << 0., 1., 1.;

    double const V1 = 0.1, V2 = 0.1, V3 = 0.1;
    double h = 2.;

    sbs::math::poly6_kernel_t W1(X1, h), W2(X2, h), W3(X3, h);
    Eigen::Vector3d gradW1 = W1.grad(Xi);
    Eigen::Vector3d gradW2 = W2.grad(Xi);
    Eigen::Vector3d gradW3 = W3.grad(Xi);

    Eigen::Matrix3d const L = V1 * gradW1 * (X1 - Xi).transpose() +
                              V2 * gradW2 * (X2 - Xi).transpose() +
                              V3 * gradW3 * (X3 - Xi).transpose();

    double const det = L.determinant();

    return 0;
}