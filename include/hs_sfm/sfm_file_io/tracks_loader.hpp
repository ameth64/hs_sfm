#ifndef _HS_SFM_SFM_FILE_IO_TRACKS_LOADER_HPP_
#define _HS_SFM_SFM_FILE_IO_TRACKS_LOADER_HPP_

#include <string>
#include <map>
#include <fstream>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{
namespace fileio
{

struct TracksLoader
{
  typedef hs::sfm::TrackContainer TrackContainer;
  typedef hs::sfm::ViewInfoIndexer ViewInfoIndexer;
  typedef hs::sfm::ViewInfo ViewInfo;
  typedef int Err;

  Err operator() (const std::string& tracks_path,
                  TrackContainer& tracks,
                  ViewInfoIndexer& view_info_indexer) const
  {
    std::ifstream tracks_file(tracks_path.c_str(), std::ios::in);
    if (!tracks_file.is_open())
    {
      return -1;
    }

    size_t number_of_tracks;
    tracks_file>>number_of_tracks;
    tracks.resize(number_of_tracks);
    std::map<std::pair<size_t, size_t>, int> flags;
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views;
      tracks_file>>number_of_views;
      tracks[i].resize(number_of_views);
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t image_id, key_id;
        int flag;
        tracks_file>>image_id>>key_id>>flag;
        tracks[i][j].first = image_id;
        tracks[i][j].second = key_id;
        flags[std::make_pair(i, image_id)] = flag;
      }
    }
    view_info_indexer.SetViewInfoByTracks(tracks);

    auto itr_flag = flags.begin();
    auto itr_flag_end = flags.end();
    for (; itr_flag != itr_flag_end; ++itr_flag)
    {
      size_t track_id = itr_flag->first.first;
      size_t image_id = itr_flag->first.second;
      int flag = itr_flag->second;
      ViewInfo* view_info =
        view_info_indexer.GetViewInfoByTrackImage(track_id, image_id);
      view_info->is_blunder = (flag == 1);
    }

    return 0;
  }
};

}
}
}

#endif
