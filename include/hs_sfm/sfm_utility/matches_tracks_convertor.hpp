#ifndef _HS_SFM_SFM_UTILITY_MATCHES_TRACKS_CONVERTOR_HPP_
#define _HS_SFM_SFM_UTILITY_MATCHES_TRACKS_CONVERTOR_HPP_

#include <map>

#include "hs_sfm/sfm_utility/match_type.hpp"

namespace hs
{
namespace sfm
{

struct MatchesTracksConvertor
{
  typedef int Err;

  Err operator() (const MatchContainer& matches,
                  TrackContainer& tracks) const
  {
    typedef std::pair<size_t, size_t> View;
    typedef std::map<View, size_t> ViewTrackIndexer;
    tracks.clear();
    MatchContainer::const_iterator itr_image_pair = matches.begin();
    MatchContainer::const_iterator itr_image_pair_end = matches.end();
    ViewTrackIndexer view_track_indexer;
    size_t number_of_tracks = 0;
    for (; itr_image_pair != itr_image_pair_end; ++itr_image_pair)
    {
      size_t image0_id = itr_image_pair->first.first;
      size_t image1_id = itr_image_pair->first.second;
      KeyPairContainer::const_iterator itr_key_pair =
        itr_image_pair->second.begin();
      KeyPairContainer::const_iterator itr_key_pair_end =
        itr_image_pair->second.end();
      for (; itr_key_pair != itr_key_pair_end; ++itr_key_pair)
      {
        size_t key0_id = itr_key_pair->first;
        size_t key1_id = itr_key_pair->second;
        View view0(image0_id, key0_id);
        View view1(image1_id, key1_id);
        ViewTrackIndexer::iterator itr_view0 = view_track_indexer.find(view0);
        ViewTrackIndexer::iterator itr_view1 = view_track_indexer.find(view1);
        if (itr_view0 != view_track_indexer.end())
        {
          if (itr_view1 != view_track_indexer.end())
          {
            if (itr_view0->second != itr_view1->second)
            {
              return -1;
            }
          }
          else
          {
            size_t track_id = itr_view0->second;
            itr_view1->second = track_id;
            tracks[track_id].push_back(view1);
          }
        }
        else
        {
          if (itr_view1 != view_track_indexer.end())
          {
            size_t track_id = itr_view1->second;
            itr_view0->second = track_id;
            tracks[track_id].push_back(view0);
          }
          else
          {
            itr_view0->second = number_of_tracks;
            itr_view1->second = number_of_tracks;
            number_of_tracks++;
            Track track;
            track.push_back(view0);
            track.push_back(view1);
            tracks.push_back(track);
          }
        }
      }
    }

    return 0;
  }

  Err operator() (const TrackContainer& tracks,
                  MatchContainer& matches) const
  {
    matches.clear();
    size_t number_of_tracks = tracks.size();
    for (size_t i = 0; i < number_of_tracks; i++)
    {
      size_t number_of_views = tracks[i].size();
      for (size_t j = 0; j < number_of_views; j++)
      {
        size_t image0_id = tracks[i][j].first;
        size_t key0_id = tracks[i][j].second;
        for (size_t k = j + 1; k < number_of_views; k++)
        {
          size_t image1_id = tracks[i][k].first;
          size_t key1_id = tracks[i][k].second;
          matches[std::make_pair(image0_id, image1_id)].push_back(
            std::make_pair(key0_id, key1_id));
        }
      }
    }
    return 0;
  }
};

}
}

#endif
