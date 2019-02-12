#!/usr/bin/env python
"""
Slack interface.
"""
__author__ = "Alex Drlica-Wagner"

import os
import logging
from datetime import datetime
import pandas as pd

from obztak.utils.database import Database

def slack_post(text, token=None, channel=None):
    """Post text to slack."""
    from slackclient import SlackClient
    slack_token   = os.environ["SLACK_API_TOKEN"] if token is None else token
    slack_channel = os.environ["SLACK_API_CHANNEL"] if channel is None else channel

    sc = SlackClient(slack_token)
    kwargs = dict(method="chat.postMessage",channel=slack_channel,text=text)
    ret = sc.api_call(**kwargs)
    return ret

def slack_qcinv(token=None, channel=None, propid=None, timedelta=None, debug=False):
    """Post inventory results to Slack.
    
    Parameters:
    -----------
    token  : slack bot token
    channel: slack channel
    propid : proposal id
    timedelta: time to collect exposure info
    debug  : execute but don't post

    Returns:
    --------
    df,pkg : data frame and text package
    """
    db = Database()
    db.connect()
    df = db.qcInv(propid=propid,timedelta=timedelta)

    if not len(df):
        logging.debug("No exposures found")
        return df,None

    kwargs = dict(index=False, float_format='{:.2f}'.format, justify='right')
    package = """Observing update from `delve-bot` @ {time} CST:
```
{content} 
```""".format(content=df.fillna('').to_string(**kwargs),
              time=datetime.now().strftime("%H:%M")
              )

    logging.debug(package)

    if debug:
        logging.debug("Exiting without posting.")
        return df,package

    ret = slack_post(package, token=token, channel=channel)
    return df,package
