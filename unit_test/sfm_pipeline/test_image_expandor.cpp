#include <fstream>

#include <gtest/gtest.h>
#include <cereal/types/map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/archives/portable_binary.hpp>

#include "hs_sfm/sfm_file_io/scene_ply_saver.hpp"
#include "hs_sfm/sfm_pipeline/image_expandor.hpp"

namespace
{

TEST(TestImageExpandor, FaileCaseTest)
{
  typedef double Scalar;
  typedef hs::sfm::pipeline::ImageExpandor<Scalar> ImageExpandor;
  typedef ImageExpandor::ImageKeysetContainer ImageKeysetContainer;
  typedef ImageExpandor::IntrinsicParamsContainer IntrinsicParamsContainer;
  typedef ImageExpandor::ImageIntrinsicMap ImageIntrinsicMap;
  typedef ImageExpandor::PointContainer PointContainer;
  typedef ImageExpandor::TrackPointMap TrackPointMap;
  typedef ImageExpandor::ImageViewTracksContainer ImageViewTracksContainer;
  typedef ImageExpandor::ImageExtrinsicMap ImageExtrinsicMap;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef ImageExpandor::ExtrinsicParamsContainer ExtrinsicParamsContainer;
  typedef hs::sfm::fileio::ScenePLYSaver<Scalar, size_t> ScenePLYSaver;
  typedef ScenePLYSaver::Image Image;
  typedef ScenePLYSaver::ImageContainer ImageContainer;

  ImageExpandor image_expandor(16, Scalar(16));

  ImageKeysetContainer image_keysets;
  IntrinsicParamsContainer intrinsic_params_set;
  ImageIntrinsicMap image_intrinsic_map;
  PointContainer points;
  TrackPointMap track_point_map;
  ImageViewTracksContainer image_view_tracks_set;
  ImageExtrinsicMap image_extrinsic_map;
  ViewInfoIndexer view_info_indexer;
  ExtrinsicParamsContainer new_extrinsic_params_set;
  std::vector<size_t> new_image_ids;

  {
    std::ifstream image_keysets_file("image_keysets_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(image_keysets_file);
    archive(image_keysets);
  }
  {
    std::ifstream intrinsic_params_set_file("intrinsic_params_set_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(intrinsic_params_set_file);
    archive(intrinsic_params_set);
  }
  {
    std::ifstream image_intrinsic_map_file("image_intrinsic_map_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(image_intrinsic_map_file);
    archive(image_intrinsic_map);
  }
  {
    std::ifstream points_file("points_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(points_file);
    archive(points);
  }
  {
    std::ifstream track_point_map_file("track_point_map_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(track_point_map_file);
    archive(track_point_map);
  }
  {
    std::ifstream image_view_tracks_set_file("image_view_tracks_set_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(image_view_tracks_set_file);
    archive(image_view_tracks_set);
  }
  {
    std::ifstream image_extrinsic_map_file("image_extrinsic_map_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(image_extrinsic_map_file);
    archive(image_extrinsic_map);
  }
  {
    std::ifstream view_info_indexer_file("view_info_indexer_227.bin", std::ios::binary);
    cereal::PortableBinaryInputArchive archive(view_info_indexer_file);
    archive(view_info_indexer);
  }

  int result = image_expandor(image_keysets,
                              intrinsic_params_set,
                              image_intrinsic_map,
                              points,
                              track_point_map,
                              image_view_tracks_set,
                              image_extrinsic_map,
                              view_info_indexer,
                              new_extrinsic_params_set,
                               new_image_ids);
  ASSERT_EQ(0, result);

  IntrinsicParamsContainer intrinsic_params_set_output;
  ExtrinsicParamsContainer extrinsic_params_set_output;
  ImageContainer images_output;;
  ScenePLYSaver saver(1);

  for (size_t i = 0; i < new_image_ids.size(); i++)
  {
    size_t image_id = new_image_ids[i];
    intrinsic_params_set_output.push_back(
      intrinsic_params_set[image_intrinsic_map[image_id]]);
    extrinsic_params_set_output.push_back(
      new_extrinsic_params_set[i]);
    Image image;
    image.m_width = 6000;
    image.m_height = 4000;
    images_output.push_back(image);
  }
  saver("image_expandor_scene.ply",
        intrinsic_params_set_output,
        extrinsic_params_set_output,
        images_output,
        points);

}

}