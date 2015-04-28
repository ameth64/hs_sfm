#ifndef _HS_SFM_UNIT_TEST_TRIANGULATE_TEST_MULTIPLE_VIEW_MAXIMUM_ESTIMATOR_REAL_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_TRIANGULATE_TEST_MULTIPLE_VIEW_MAXIMUM_ESTIMATOR_REAL_DATA_GENERATOR_HPP_

#include <string>
#include <fstream>
#include <sstream>

#include "hs_math/random/normal_random_var.hpp"
#include "hs_math/random/uniform_random_var.hpp"

#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"
#include "hs_sfm/sfm_utility/projective_functions.hpp"
#include "hs_sfm/sfm_utility/matches_tracks_convertor.hpp"

namespace hs
{
namespace sfm
{
namespace triangulate
{

template <typename _Scalar>
class TestMultipleViewMaximumLikelihoodEstimatorRealDataGenerator
{
public:
  typedef _Scalar Scalar;
  typedef size_t ImageDimension;
  typedef int Err;
  typedef ImageKeys<Scalar> Keyset;
  typedef EIGEN_STD_VECTOR(Keyset) KeysetContainer;
  typedef typename Keyset::Key Key;
  typedef EIGEN_VECTOR(Scalar, 3) Point;
  typedef EIGEN_STD_VECTOR(Point) PointContainer;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef hs::sfm::ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;

private:
  typedef EIGEN_MATRIX(Scalar, 3, 3) Matrix33;

public:
  TestMultipleViewMaximumLikelihoodEstimatorRealDataGenerator(
    Scalar outlier_ratio, Scalar key_stddev)
    : outlier_ratio_(outlier_ratio), key_stddev_(key_stddev) {}

  Err Generate(const std::string& bundler_out_path,
               const IntrinsicParamsContainer& intrinsic_params_set,
               const hs::sfm::ObjectIndexMap& image_intrinsic_map,
               ExtrinsicParamsContainer& extrinsic_params_set,
               hs::sfm::ObjectIndexMap& image_extrinsic_map,
               ImageDimension image_width, ImageDimension image_height,
               hs::sfm::TrackContainer& tracks,
               KeysetContainer& keysets_noised,
               PointContainer& points_true) const
  {
    Err result = 0;
    while (1)
    {
      KeysetContainer keysets_true;
      result = LoadBundlerOutFile(bundler_out_path,
                                  intrinsic_params_set,
                                  image_intrinsic_map,
                                  image_extrinsic_map,
                                  extrinsic_params_set,
                                  points_true,
                                  keysets_true,
                                  tracks);
      if (result != 0) break;

      Image image;
      image.m_width = image_width;
      image.m_height = image_height;
      ImageContainer images(extrinsic_params_set.size(), image);

      result = GenerateNoisedKeysets(keysets_true,
                                     images,
                                     keysets_noised);
      if (result != 0) break;
      break;
    }
    return result;
  }

private:
  Err LoadBundlerOutFile(
    const std::string& bundler_out_path,
    const IntrinsicParamsContainer& intrinsic_params_set,
    const hs::sfm::ObjectIndexMap& image_intrinsic_map,
    hs::sfm::ObjectIndexMap& image_extrinsic_map,
    ExtrinsicParamsContainer& extrinsic_params_set,
    PointContainer& points_true,
    KeysetContainer& keysets_true,
    hs::sfm::TrackContainer& tracks) const
  {
    std::ifstream bundler_out_file(bundler_out_path);
    if (!bundler_out_file)
    {
      return -1;
    }

    std::stringstream ss;
    std::string line;
    std::getline(bundler_out_file, line);
    std::getline(bundler_out_file, line);
    ss.str(line);
    size_t number_of_images;
    size_t number_of_points;
    ss>>number_of_images>>number_of_points;

#ifdef GETLINE
#undef GETLINE
#endif
#define GETLINE \
        std::getline(bundler_out_file, line);\
        ss.clear();\
        ss.str(line);

    image_extrinsic_map.Resize(number_of_images);
    for (size_t image_id = 0; image_id < number_of_images; image_id++)
    {
      GETLINE;
      Scalar focal;
      ss>>focal;
      if (focal == 0) continue;
      Matrix33 rotation_matrix;
      GETLINE;
      ss>>rotation_matrix(0, 0)>>rotation_matrix(0, 1)>>rotation_matrix(0, 2);
      GETLINE;
      ss>>rotation_matrix(1, 0)>>rotation_matrix(1, 1)>>rotation_matrix(1, 2);
      GETLINE;
      ss>>rotation_matrix(2, 0)>>rotation_matrix(2, 1)>>rotation_matrix(2, 2);
      GETLINE;
      Point translate;
      ss>>translate[0]>>translate[1]>>translate[2];
      ExtrinsicParams extrinsic_params;
      extrinsic_params.position() = -rotation_matrix.transpose() * translate;
      rotation_matrix.row(1) *= Scalar(-1);
      rotation_matrix.row(2) *= Scalar(-1);
      extrinsic_params.rotation() = rotation_matrix;
      extrinsic_params_set.push_back(extrinsic_params);
      image_extrinsic_map[image_id] = extrinsic_params_set.size() - 1;
    }

    keysets_true.resize(extrinsic_params_set.size());
    for (size_t point_id = 0; point_id < number_of_points; point_id++)
    {
      GETLINE;
      Point point;
      ss>>point[0]>>point[1]>>point[2];
      GETLINE;
      GETLINE;
      size_t number_of_views;
      ss>>number_of_views;
      hs::sfm::Track track;
      for (size_t view_id = 0; view_id < number_of_views; view_id++)
      {
        size_t image_id, key_id;
        Scalar key_x, key_y;
        ss>>image_id>>key_id>>key_x>>key_y;
        key_x = 3000 + key_x;
        key_y = 2000 - key_y;
        if (image_extrinsic_map.IsValid(image_id))
        {
          size_t extrinsic_id = image_extrinsic_map[image_id];
          const ExtrinsicParams& extrinsic_params =
            extrinsic_params_set[extrinsic_id];
          size_t intrinsic_id = image_intrinsic_map[image_id];
          const IntrinsicParams& intrinsic_params =
            intrinsic_params_set[intrinsic_id];

          Key key = ProjectiveFunctions<Scalar>::WorldPointProjectToImageKey(
              intrinsic_params, extrinsic_params, point);

          keysets_true[extrinsic_id].AddKey(key);
          track.push_back(
            std::make_pair(image_id, keysets_true[extrinsic_id].size() - 1));
        }
      }

      if (track.size() > 1)
      {
        tracks.push_back(track);
        points_true.push_back(point);
      }
    }

    //Point point_mean = Point::Zero();
    //for (size_t point_id = 0; point_id < points_true.size(); point_id++)
    //{
    //  point_mean += points_true[point_id];
    //}
    //point_mean /= Scalar(points_true.size());
    //for (size_t point_id = 0; point_id < points_true.size(); point_id++)
    //{
    //  points_true[point_id] -= point_mean;
    //}
    //for (size_t extrinsic_id = 0; extrinsic_id < extrinsic_params_set.size();
    //     extrinsic_id++)
    //{
    //  extrinsic_params_set[extrinsic_id].position() -= point_mean;
    //}

#undef GETLINE

    return 0;
  }

  Err GenerateNoisedKeysets(const KeysetContainer& keysets_true,
                            const ImageContainer& images,
                            KeysetContainer& keysets_noised,
                            bool generate_outliers = true) const
  {
    keysets_noised = keysets_true;
    EIGEN_MATRIX(Scalar, 2, 2) covariance_key;
    covariance_key.setIdentity();
    covariance_key *= key_stddev_ * key_stddev_;
    size_t number_of_keysets = keysets_true.size();
    for (size_t i = 0; i < number_of_keysets; i++)
    {
      size_t number_of_keys = keysets_true[i].size();
      for (size_t j = 0; j < number_of_keys; j++)
      {
        Scalar random;
        hs::math::random::UniformRandomVar<Scalar, 1>::Generate(
          0, 1, random);
        if (random < outlier_ratio_ && generate_outliers)
        {
          Key min_key;
          min_key.setZero();
          Key max_key;
          max_key << Scalar(images[i].m_width),
                     Scalar(images[i].m_height);
          hs::math::random::UniformRandomVar<Scalar, 2>::Generate(
            min_key, max_key, keysets_noised[i][j]);
        }
        else
        {
          Key mean = keysets_true[i][j];
          hs::math::random::NormalRandomVar<Scalar, 2>::Generate(
            mean, covariance_key, keysets_noised[i][j]);
        }
      }
    }

    return 0;
  }

private:
  Scalar outlier_ratio_;
  Scalar key_stddev_;
};

}
}
}

#endif