#ifndef _HS_SFM_UTILITY_MATCHE_TYPE_HPP_
#define _HS_SFM_UTILITY_MATCHE_TYPE_HPP_

#include <vector>
#include <map>

#include "cereal/types/vector.hpp"
#include "cereal/types/map.hpp"
#include "cereal/types/utility.hpp"

namespace hs
{
namespace sfm
{

typedef std::vector<std::vector<int> > MatchGuide;
typedef std::pair<size_t, size_t> ImagePair;
typedef std::pair<size_t, size_t> KeyPair;
typedef std::vector<KeyPair> KeyPairContainer;
typedef std::map<ImagePair, KeyPairContainer> MatchContainer;
typedef std::vector<std::pair<size_t, size_t> > Track;
typedef std::vector<Track> TrackContainer;
typedef std::vector<std::pair<size_t, size_t> > CameraView;
typedef std::vector<CameraView> CameraViewContainer;

class ObjectIndexMap	//记录每张图像的?
{
public:
  ObjectIndexMap(){}
  ObjectIndexMap(const ObjectIndexMap& other) {mapper_ = other.mapper_;}
  ObjectIndexMap(size_t number_of_objects)
    : mapper_(number_of_objects, invalid_value()) {}
public:

  ObjectIndexMap& operator = (const ObjectIndexMap& other)
  {
    mapper_ = other.mapper_;
    return *this;
  }

  size_t GetMappedId(size_t object_id) const
  {
    return mapper_[object_id];
  }

  void SetObjectId(size_t object_id, size_t mapped_id)
  {
    mapper_[object_id] = mapped_id;
  }

  size_t operator[] (size_t object_id) const
  {
    return mapper_[object_id];
  }

  size_t& operator[] (size_t object_id)
  {
    return mapper_[object_id];
  }

  bool IsValid(size_t object_id) const
  {
    return (mapper_[object_id] != invalid_value());
  }

  void Resize(size_t number_of_objects)
  {
    mapper_.resize(number_of_objects, invalid_value());
  }

  size_t Size() const
  {
    return mapper_.size();
  }

  void AddObject(size_t mapped_id)
  {
    mapper_.push_back(mapped_id);
  }

  bool operator == (const ObjectIndexMap& other) const
  {
    return (mapper_ == other.mapper_);
  }

  template <typename Archive>
  void save(Archive& archive) const
  {
    archive(mapper_);
  }

  template <typename Archive>
  void load(Archive& archive)
  {
    archive(mapper_);
  }

private:
  inline size_t invalid_value() const
  {
    return std::numeric_limits<size_t>::max();
  }

private:
  std::vector<size_t> mapper_;
};

struct ViewInfo
{
  size_t track_id;
  size_t image_id;
  size_t key_id;
  bool is_blunder;

  bool operator == (const ViewInfo& other) const
  {
    return (track_id == other.track_id &&
      image_id == other.image_id &&
      key_id == other.key_id &&
      is_blunder == other.is_blunder);
  }

  template <typename Archive>
  void serialize(Archive& archive)
  {
    archive(track_id, image_id, key_id, is_blunder);
  }
};

class ViewInfoIndexer
{
public:
  void Clear()
  {
    views_info_.clear();
    track_image_index_.clear();
    image_key_index_.clear();
  }

  void SetViewInfoByTracks(const hs::sfm::TrackContainer& tracks)
  {
    Clear();
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      for (size_t j = 0; j < number_of_views; j++)
      {
        ViewInfo view_info;
        view_info.track_id = i;
        view_info.image_id = tracks[i][j].first;
        view_info.key_id = tracks[i][j].second;
        view_info.is_blunder = false;
        views_info_.push_back(view_info);
        size_t view_id = views_info_.size() - 1;
        track_image_index_[std::make_pair(view_info.track_id,
                                          view_info.image_id)] = view_id;
        image_key_index_[std::make_pair(view_info.image_id,
                                        view_info.key_id)] = view_id;
      }
    }
  }

  const ViewInfo* GetViewInfoByTrackImage(size_t track_id,
                                          size_t image_id) const
  {
    auto view_itr = track_image_index_.find(std::make_pair(track_id,
                                                            image_id));
    if (view_itr != track_image_index_.end())
    {
      return &views_info_[view_itr->second];
    }
    else
    {
      return nullptr;
    }
  }

  ViewInfo* GetViewInfoByTrackImage(size_t track_id,
                                    size_t image_id)
  {
    auto view_itr = track_image_index_.find(std::make_pair(track_id,
                                                            image_id));
    if (view_itr != track_image_index_.end())
    {
      return &views_info_[view_itr->second];
    }
    else
    {
      return nullptr;
    }
  }

  const ViewInfo* GetViewInfoByImageKey(size_t image_id,
                                           size_t key_id) const
  {
    auto view_itr = image_key_index_.find(std::make_pair(image_id, key_id));
    if (view_itr == image_key_index_.end())
    {
      return nullptr;
    }
    else
    {
      return &views_info_[view_itr->second];
    }
  }

  ViewInfo* GetViewInfoByImageKey(size_t image_id,
                                     size_t key_id)
  {
    auto view_itr = image_key_index_.find(std::make_pair(image_id, key_id));
    if (view_itr == image_key_index_.end())
    {
      return nullptr;
    }
    else
    {
      return &views_info_[view_itr->second];
    }
  }

  bool operator == (const ViewInfoIndexer& other) const
  {
    return (views_info_ == other.views_info_ &&
            track_image_index_ == other.track_image_index_ &&
            image_key_index_ == other.image_key_index_);
  }

  template <typename Archive>
  void serialize(Archive& archive)
  {
    archive(views_info_, track_image_index_, image_key_index_);
  }

private:
  std::vector<ViewInfo> views_info_;
  std::map<std::pair<size_t, size_t>, size_t> track_image_index_;
  std::map<std::pair<size_t, size_t>, size_t> image_key_index_;
};


}
}

#endif
