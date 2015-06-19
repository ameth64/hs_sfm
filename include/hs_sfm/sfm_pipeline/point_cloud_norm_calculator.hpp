#ifndef _HS_SFM_SFM_PIPELINE_POINT_CLOUD_NORM_CALCULATOR_HPP_
#define _HS_SFM_SFM_PIPELINE_POINT_CLOUD_NORM_CALCULATOR_HPP_

#include <utility>
#include <limits>

#include <flann/flann.hpp>

#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class PointCloudNormCalculator
{
public:
  typedef _Scalar Scalar;
  typedef CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_VECTOR(Scalar, 3) Point;

  struct CameraParams
  {
    IntrinsicParams intrinsic_params;
    ExtrinsicParams extrinsic_params;
    size_t image_width;
    size_t image_height;
  };

  typedef EIGEN_STD_VECTOR(CameraParams) CameraParamsContainer;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
private:
  typedef flann::L2<Scalar> Distance;
  typedef flann::Index<Distance> Index;
  typedef EIGEN_VECTOR(Scalar, 2) Key;

public:
  PointCloudNormCalculator(int knn = 64)
    : knn_(knn) {}

  int operator() (const CameraParamsContainer& cameras,
                  const PointContainer& points,
                  PointContainer& norms) const
  {
    Index index = BuildQueryTree(points);
    norms.resize(points.size());
    for (size_t i = 0; i < points.size(); i++)
    {
      Point norm = PCANearestNeighborNorm(points, i, index);

      //计算所有与该点关联的相机中，相机拍摄方向与法向量的夹角
      Scalar max_angle = -std::numeric_limits<Scalar>::max();
      Scalar max_angle_signed = Scalar(0);
      for (size_t j = 0; j < cameras.size(); j++)
      {
        Key key =
          ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
            cameras[j].intrinsic_params,
            cameras[j].extrinsic_params,
            points[i]);
        if (key[0] >= Scalar(0) && key[0] < Scalar(cameras[j].image_width) &&
            key[1] >= Scalar(0) && key[1] < Scalar(cameras[j].image_height))
        {
          EIGEN_MATRIX(Scalar, 3, 3) r_matrix =
            cameras[j].extrinsic_params.rotation();
          Point camera_direction = r_matrix.row(2).transpose();
          Scalar angle = camera_direction.dot(norm);
          if (max_angle < std::abs(angle))
          {
            max_angle = std::abs(angle);
            max_angle_signed = angle;
          }
        }
      }

      if (max_angle_signed > Scalar(0))
      {
        norm *= Scalar(-1);
      }

      norms[i] = norm;
    }
    return 0;
  }

  int operator() (const PointContainer& points,
    const Point& up_vector,
    PointContainer& norms) const
  {
    Index index = BuildQueryTree(points);
    norms.resize(points.size());
    for (size_t i = 0; i < points.size(); i++)
    {
      Point norm = PCANearestNeighborNorm(points, i, index);

      if (norm.dot(up_vector) < Scalar(0))
      {
        norm *= Scalar(-1);
      }

      norms[i] = norm;
    }
    return 0;
  }

private:
  Index BuildQueryTree(const PointContainer& points) const
  {
    std::vector<Scalar> data(points.size() * 3);
    for (size_t i = 0; i < points.size(); i++)
    {
      data[i * 3 + 0] = points[i][0];
      data[i * 3 + 1] = points[i][1];
      data[i * 3 + 2] = points[i][2];
    }

    flann::Matrix<Scalar> index_matrix(data.data(), points.size(), 3);
    Index index(index_matrix, flann::KDTreeSingleIndexParams());
    index.buildIndex();
    return index;
  }

  Point PCANearestNeighborNorm(const PointContainer& points,
                               size_t id,
                               const Index& index) const
  {
    typedef typename Distance::ResultType DistanceResultType;
    typedef EIGEN_MATRIX(Scalar, 3, Eigen::Dynamic) Matrix3X;
    typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

    Point norm;
    Scalar search_data[3] = { points[id][0], points[id][1], points[id][2] };
    std::vector<size_t> indices_data(knn_);
    std::vector<DistanceResultType> distances_data(knn_);

    flann::Matrix<Scalar> serach_matrix(search_data, 1, 3);
    flann::Matrix<size_t> indices_matrix(indices_data.data(), 1, knn_);
    flann::Matrix<DistanceResultType> distances_matrix(
      distances_data.data(), 1, knn_);
    index.knnSearch(serach_matrix, indices_matrix, distances_matrix, knn_,
                    flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

    Point mean = Point::Zero();
    for (size_t j = 0; j < knn_; j++)
    {
      mean[0] += points[indices_matrix[0][j]][0];
      mean[1] += points[indices_matrix[0][j]][1];
      mean[2] += points[indices_matrix[0][j]][2];
    }
    mean /= Scalar(knn_);

    Matrix3X A(3, knn_);
    for (size_t j = 0; j < knn_; j++)
    {
      Point point_neighbor = points[indices_matrix[0][j]];
      A.col(j) = point_neighbor - mean;
    }
    Matrix33 cov = A * A.transpose();
    Eigen::EigenSolver<Matrix33> eigen_solver(cov);
    Scalar min_eigen_value = std::numeric_limits<Scalar>::max();
    int min_idx = -1;
    for (int j = 0; j < 3; j++)
    {
      Scalar eigen_value = eigen_solver.eigenvalues()[j].real();
      if (min_eigen_value > eigen_value)
      {
        min_eigen_value = eigen_value;
        min_idx = j;
      }
    }

    norm = eigen_solver.eigenvectors().col(min_idx).real();
    norm /= norm.norm();

    return norm;
  }

private:
  int knn_;
};

}
}
}

#endif