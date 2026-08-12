"""Constants for the Global Meteor Network observation layer."""

from datetime import timedelta

GMN_DAILY_INDEX_URL = (
    "https://globalmeteornetwork.org/data/traj_summary_data/daily/"
)
GMN_CACHE_SCHEMA = "zstarview.gmn-meteor-cache.v1"
GMN_INDEX_FRESH_TTL = timedelta(hours=6)
GMN_RECENT_FILE_FRESH_TTL = timedelta(hours=6)
GMN_WINDOW = timedelta(hours=24)
GMN_LATEST_LOOKBACK = timedelta(days=7)
GMN_CANDIDATE_RADIUS_KM = 500.0
GMN_MAX_DISPLAY_TRAILS = 100
