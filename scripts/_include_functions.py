#!/usr/bin/env python
# 🐛

# Function to send myself push alerts when things take a long time (python version)

import pushover
import platform
import os
import pwd

# User key and token must be set (omitted below)
ukey = ''
token = ''

def push_status(topic=None,msg=None,user_key=ukey,token=token):
	if topic is None:
		topic = 'Action'
		subj = 'python-pushover'
	else:
		subj = 'python-pushover - '+topic
	if msg is None:
		m_msg = pwd.getpwuid(os.getuid())[0]+'@'+platform.node()+': '+topic+' completed'
	else:
		m_msg = msg
	pushover.Client(user_key,api_token=token).send_message(m_msg,title=subj)
