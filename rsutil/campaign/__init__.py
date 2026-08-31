"""Deterministic, lightweight campaign specification and generation.

Campaigns describe setup inputs and external calculation metadata.  They do
not execute RS-LMTO jobs or provide workflow/scheduler functionality.
"""

from .generation import (
    CAMPAIGN_INDEX_FILENAME,
    RUN_MANIFEST_FILENAME,
    CampaignCollisionError,
    CampaignError,
    CampaignGenerator,
    CampaignValidation,
    generate_campaign,
    load_campaign,
    validate_campaign,
)
from .models import (
    CAMPAIGN_SCHEMA_VERSION,
    Campaign,
    CampaignEntry,
    CampaignIndex,
    CampaignIndexEntry,
)

__all__ = [
    "CAMPAIGN_INDEX_FILENAME",
    "CAMPAIGN_SCHEMA_VERSION",
    "RUN_MANIFEST_FILENAME",
    "Campaign",
    "CampaignCollisionError",
    "CampaignEntry",
    "CampaignError",
    "CampaignGenerator",
    "CampaignIndex",
    "CampaignIndexEntry",
    "CampaignValidation",
    "generate_campaign",
    "load_campaign",
    "validate_campaign",
]
