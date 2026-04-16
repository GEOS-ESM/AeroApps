#!/usr/bin/env python3
"""
    Script to run a  buddy check on AMoN sites 
"""
import os, sys
import argparse
from pyobs import buddycheck
from load_amon import LOAD_AMON


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("config",help='configuration yaml file')

    args = parser.parse_args()

    df_amon, sites = LOAD_AMON(args.config,model_name='baseline')

    # rename columns to names expected by buddycheck
    df_amon = df_amon.rename(columns={
            'SITEID': 'station_id',
            'NH3': 'obs',
            'baseline': 'model',
            })

    sites = sites.rename(columns={
            'siteId': 'station_id',
            'latitude': 'lat',
            'longitude': 'lon',
            })

    # configuration parameters for the buddy check
    config = buddycheck.InnovationBuddyConfig(
                radius_km=150.0,
                min_neighbors=3,
                distance_scale_km=75.0,
                buddy_method="weighted_mean",
                outlier_sigma_thresh=2.5,
                min_valid_buddy_count=10,
    )

    buddy_time_df, station_scores_df = buddycheck.score_stations_innovation_based(
        obs_df=df_amon,
        station_meta=sites,
        config=config,
    )

    # implement an alternative review flag
    review_flag = ((station_scores_df["suspect_score"] >= 6.0) &
                  (
                    (station_scores_df["z_f_bias_vs_buddies"] >= 2.5) |
                    (station_scores_df["z_f_rmse_vs_buddies"] >= 2.5) |
                    (station_scores_df["z_f_frac_outlier"] >= 2.5)
                  )
                 )

    station_scores_df["review_flag_new"] = review_flag

    # print the top 20 stations
    print(station_scores_df.sort_values("suspect_score", ascending=False).head(20))

    files = buddycheck.plot_flagged_stations(
        station_scores_df=station_scores_df,
        buddy_time_df=buddy_time_df,
        station_meta=sites,
        outdir="plots/amon_flagged_stations",
        only_flagged=True,
#        top_n=5,
        )
