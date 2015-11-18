#ifndef _HS_SFM_SFM_FILE_IO_SCENE_PLY_SAVER_HPP_
#define _HS_SFM_SFM_FILE_IO_SCENE_PLY_SAVER_HPP_

#include <string>
#include <fstream>

#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

/**
 *  给定相机内外参数以及相片的宽度与高度，及点云数据，输出成ply格式文件。
 *
 *  该类主要方便内部调试使用，可直观的看到场景的大致状态。
 */
template <typename _Scalar, typename _ImageDimension>
class ScenePLYSaver
{
public:
  typedef _Scalar Scalar;
  typedef _ImageDimension ImageDimension;
  typedef int Err;
  typedef hs::sfm::CameraIntrinsicParams<Scalar> IntrinsicParams;
  typedef EIGEN_STD_VECTOR(IntrinsicParams) IntrinsicParamsContainer;
  typedef hs::sfm::CameraExtrinsicParams<Scalar> ExtrinsicParams;
  typedef EIGEN_STD_VECTOR(ExtrinsicParams) ExtrinsicParamsContainer;
  typedef hs::sfm::ImageParams<ImageDimension> Image;
  typedef EIGEN_STD_VECTOR(Image) ImageContainer;
  typedef EIGEN_VECTOR(Scalar, 3) Point3D;
  typedef EIGEN_STD_VECTOR(Point3D) Point3DContainer;

private:
  typedef typename ExtrinsicParams::Position Position;

public:
  ScenePLYSaver(Scalar camera_size) : camera_size_(camera_size) {}

  Err operator() (const std::string& ply_path,
                  const IntrinsicParamsContainer& intrinsic_params_set,
                  const ExtrinsicParamsContainer& extrinsic_params_set,
                  const ImageContainer& images,
                  const Point3DContainer& points,
                  const Point3DContainer* norms = nullptr) const
  {
    size_t number_of_images = intrinsic_params_set.size();
    if (number_of_images != extrinsic_params_set.size()) return -1;
    if (number_of_images != images.size()) return -1;
    size_t number_of_points = points.size();
    std::ofstream ply_file(ply_path.c_str());
    if (!ply_file) return -1;

    size_t number_of_vertices = number_of_images * 5 + number_of_points;
    ply_file << "ply\n"
             << "format ascii 1.0\n"
             << "element vertex "<<number_of_vertices<<"\n"
             << "property float x\n"
             << "property float y\n"
             << "property float z\n";

    if (norms)
    {
      ply_file << "property float nx\n"
               << "property float ny\n"
               << "property float nz\n";
    }
    else
    {
      ply_file << "element face "<<number_of_images * 5<<"\n"
               << "property list uchar int vertex_index\n";
    }
    ply_file << "end_header\n";

    for (size_t i = 0; i < number_of_images; i++)
    {
      const IntrinsicParams& intrinsic_params = intrinsic_params_set[i];
      const ExtrinsicParams& extrinsic_params = extrinsic_params_set[i];
      const Image& image = images[i];
      const Position& position = extrinsic_params.position();
      ply_file << position[0] << " " << position[1] << " " << position[2];
      if (norms)
      {
        ply_file << " 0 0 1";
      }
      ply_file << "\n";
      Position corner[4];
      Scalar width_ratio =
        Scalar(image.m_width) / intrinsic_params.focal_length() *
        camera_size_ / Scalar(2);
      Scalar height_ratio =
        Scalar(image.m_height) / intrinsic_params.focal_length() *
        camera_size_ / Scalar(2);
      corner[0] << -width_ratio,
                   -height_ratio,
                   camera_size_;
      corner[1] << width_ratio,
                   -height_ratio,
                   camera_size_;
      corner[2] << width_ratio,
                   height_ratio,
                   camera_size_;
      corner[3] << -width_ratio,
                   height_ratio,
                   camera_size_;
      EIGEN_MATRIX(Scalar, 3, 3) rotation_matrix = extrinsic_params.rotation();
      Position translate = -rotation_matrix * position;
      for (int j = 0; j < 4; j++)
      {
        corner[j] = rotation_matrix.transpose() * (corner[j] - translate);
        ply_file<<corner[j][0]<<" "<<corner[j][1]<<" "<<corner[j][2];
        if (norms)
        {
          ply_file << " 0 0 1";
        }
        ply_file << "\n";
      }
    }

    for (size_t i = 0; i < number_of_points; i++)
    {
      ply_file<<points[i][0]<<" "<<points[i][1]<<" "<<points[i][2];
      if (norms)
      {
        ply_file << " " << (*norms)[i][0]
                 << " " << (*norms)[i][1]
                 << " " << (*norms)[i][2];
      }
      ply_file << "\n";
    }

    if (!norms)
    {
      for (size_t i = 0; i < number_of_images; i++)
      {
        ply_file<<"4 "
                <<i * 5 + 1<<" "
                <<i * 5 + 2<<" "
                <<i * 5 + 3<<" "
                <<i * 5 + 4<<"\n";
        ply_file<<"3 "<<i * 5 + 0<<" "<<i * 5 + 1<<" "<<i * 5 + 2<<"\n";
        ply_file<<"3 "<<i * 5 + 0<<" "<<i * 5 + 2<<" "<<i * 5 + 3<<"\n";
        ply_file<<"3 "<<i * 5 + 0<<" "<<i * 5 + 3<<" "<<i * 5 + 4<<"\n";
        ply_file<<"3 "<<i * 5 + 0<<" "<<i * 5 + 4<<" "<<i * 5 + 1<<"\n";
      }
    }

    return 0;
  }

private:
  Scalar camera_size_;

};

}
}
}

#endif
