#ifndef _HS_SFM_UNIT_TEST_SFM_PIPELINE_SYNTHETIC_1DSFM_DATA_GENERATOR_HPP_
#define _HS_SFM_UNIT_TEST_SFM_PIPELINE_SYNTHETIC_1DSFM_DATA_GENERATOR_HPP_

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <algorithm>

#include "hs_sfm/sfm_utility/match_type.hpp"
#include "hs_sfm/sfm_utility/camera_type.hpp"

namespace hs
{
namespace sfm
{
namespace pipeline
{

template <typename _Scalar>
class Synthetic1DSFMDataGenerator
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

  Err Generate(const std::string& dataset_dir,
               IntrinsicParamsContainer& intrinsic_params_set,
               hs::sfm::ObjectIndexMap& image_intrinsic_map,
               hs::sfm::TrackContainer& tracks,
               KeysetContainer& keysets) const
  {
    Err result = 0;

    while (1)
    {
      result = LoadCoords(dataset_dir,
                          intrinsic_params_set,
                          image_intrinsic_map,
                          keysets);
      if (result != 0) break;

      result = LoadTracks(dataset_dir,
                          tracks);
      if (result != 0) break;

      break;
    }

    return result;
  }

private:
  Err LoadCoords(const std::string& dataset_dir,
                 IntrinsicParamsContainer& intrinsic_params_set,
                 hs::sfm::ObjectIndexMap& image_intrinsic_map,
                 KeysetContainer& keysets) const
  {
    std::string coords_path = dataset_dir + "/coords.txt";
    std::ifstream coords_file(coords_path);
    if (!coords_file) return -1;

    std::stringstream ss;
    std::string line;
#ifdef GETLINE
#undef GETLINE
#endif
#define GETLINE \
        std::getline(coords_file, line);\
        ss.clear();\
        ss.str(line);

    std::vector<size_t> image_ids;
    intrinsic_params_set.clear();
    KeysetContainer keysets_tight;
    while (!coords_file.eof())
    {
      GETLINE;
      if (line.empty())
      {
        break;
      }
      std::string dummy;
      size_t image_id;
      size_t number_of_keys = 0;
      Scalar px, py, focal;
      ss>>dummy>>dummy>>image_id>>dummy>>dummy>>dummy>>dummy>>dummy>>dummy
        >>number_of_keys>>dummy>>dummy>>dummy
        >>px>>dummy>>dummy>>dummy
        >>py>>dummy>>dummy>>dummy
        >>focal;
      Keyset keyset;
      for (size_t i = 0; i < number_of_keys; i++)
      {
        GETLINE;
        Key key;
        size_t key_id;
        ss>>key_id>>key[0]>>key[1];
        keyset.AddKey(key);
      }
      image_ids.push_back(image_id);
      IntrinsicParams intrinsic_params(focal, 0, px, py);
      intrinsic_params_set.push_back(intrinsic_params);
      keysets_tight.push_back(keyset);
    }

    size_t number_of_images = *(std::max_element(image_ids.begin(),
                                                 image_ids.end())) + 1;
    image_intrinsic_map.Resize(number_of_images);
    keysets.resize(number_of_images);
    for (size_t i = 0; i < image_ids.size(); i++)
    {
      image_intrinsic_map.SetObjectId(image_ids[i], i);
      keysets[image_ids[i]] = keysets_tight[i];
    }

#undef GETLINE

    return 0;
  }

  Err LoadTracks(const std::string& dataset_dir,
                 hs::sfm::TrackContainer& tracks) const
  {
    std::string tracks_path = dataset_dir + "/tracks.txt";
    std::ifstream tracks_file(tracks_path);
    if (!tracks_file) return -1;

    size_t number_of_tracks = 0;
    tracks_file>>number_of_tracks;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = 0;
      tracks_file>>number_of_views;
      hs::sfm::Track track;
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t image_id, key_id;
        tracks_file>>image_id>>key_id;
        track.push_back(std::make_pair(image_id, key_id));
      }
      tracks.push_back(track);
    }

    return 0;
  }
};

}
}
}

#endif
