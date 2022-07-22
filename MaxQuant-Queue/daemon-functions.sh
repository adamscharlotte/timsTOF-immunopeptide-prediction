#!/bin/bash

## ################################ ##
##                                  ##
## daemon-functions.sh VERSION 0.1a ##
##                                  ##
## ################################ ##

## Distrubuted under the GPL
## http://www.gnu.org/licenses/gpl-3.0.txt

## No warranty of any kind... May run 
## off with your daughter. May explode
## in a ball of smoke and fire. Might
## work. Use at your own risk

## #!/bin/bash
##
## # Example Usage: datelogger.sh
## #	A sample daemon which simply logs the
## #	date and time once per second.
##
## function payload() {
##   while [ true ]; do
##     checkforterm
##     date
##     sleep 1
##   done
## }
##
## source /path/to/daemon-functions.sh

function daemonize() {
	echo $MY_PID > $MY_PIDFILE
	exec 3>&-           # close stdin
	exec 2>>$MY_ERRFILE # redirect stderr
	exec 1>>$MY_LOGFILE # redirect stdout
	echo $(date)" Daemonizing $MY_PID $HOSTNAME" >> $MY_ERRFILE
}

function checkforterm() {
	if [ -f $MY_KILLFILE ]; then
		echo $(date)" Terminating gracefully" >> $MY_ERRFILE
		rm $MY_PIDFILE
		rm $MY_KILLFILE
		kill $MY_PID
		exit 0
	fi
	sleepcount=0
	while [ -f $MY_WAITFILE ]; do 
		let sleepcount=$sleepcount+1
		let pos=$sleepcount%10
		if [ $pos -eq 0 ]; then
			echo $(date)" Sleeping..."
			echo $(date)" Sleeping..." >> $MY_ERRFILE
		fi
		if [ -f $MY_KILLFILE ]; then
			rm $MY_WAITFILE
			checkforterm
		fi
		sleep 1
	done
}

MY_PID=$$
MY_PATH=$(readlink -f $0)
MY_ROOT=$(dirname $MY_PATH)
MY_NAME=$(basename $MY_PATH)
MY_PIDFILE="$MY_ROOT/.$MY_NAME-$HOSTNAME.pid"
MY_KILLFILE="$MY_ROOT/.$MY_NAME-$HOSTNAME.kill"
MY_ERRFILE="$MY_ROOT/.$MY_NAME-$HOSTNAME.err"
MY_LOGFILE="$MY_ROOT/.$MY_NAME-$HOSTNAME.log"
MY_WAITFILE="$MY_ROOT/.$MY_NAME-$HOSTNAME.wait"
MY_BLOCKFILE="$MY_ROOT/.$MY_NAME-$HOSTNAME.block"

CR="
"
SP=" "
OIFS=$IFS

case $1 in
	pause)
		touch $MY_WAITFILE
		;;
	resume)
		rm $MY_WAITFILE
		;;
	restart)
		$0 stop
		$0 start
		;;
	start)
		if [ -f $MY_BLOCKFILE ]; then
			echo "Daemon execution has been disabled"
			exit 0
		fi
		$0 run &
		echo "Daemon Started"
		exec 3>&- # close stdin
		exec 2>&- # close stderr
		exec 1>&- # close stdout
		exit 0
		;;
	disable)
		touch $MY_BLOCKFILE
		$0 stop
		;;
	enable)
		if [ -f $MY_BLOCKFILE ]; then rm $MY_BLOCKFILE; fi
		;;
	stop)
		echo -n "Terminating daemon... "
		$0 stat 1>/dev/null 2>/dev/null
		if [ $? -ne 0 ]; then
			echo "process is not running"
			exit 0
		fi
		touch $MY_KILLFILE
		$0 stat 1>/dev/null 2>/dev/null
		ECODE=$?
		waitcount=0
		if [ "$waitcountmax" = "" ]; then waitcountmax=30; fi
		while [ $ECODE -eq 0 ]; do
			sleep 1
			let waitcount=$waitcount+1
			if [ $waitcount -lt $waitcountmax ]; then
				$0 stat 1>/dev/null 2>/dev/null
				ECODE=$?
			else
				ECODE=1
			fi
		done
		$0 stat 1>/dev/null 2>/dev/null
		if [ $? -eq 0 ]; then
			PID=$(cat $MY_PIDFILE)
			kill $PID
			rm $MY_PIDFILE
			rm $MY_KILLFILE
			echo "Process Killed"
			echo $(date)" Terminating forcefully" >> $MY_ERRFILE
			exit 0;
		else
			echo "Process exited gracefully"
		fi
		;;
	stat|status)
		if [ -f $MY_BLOCKFILE ]; then
			echo "Daemon execution disabled"
		fi
		if [ ! -f $MY_PIDFILE ]; then
			echo "$MY_NAME is not running"
			exit 1
		fi
		pgrep -l -f "$MY_NAME run" | grep -q -E "^$(cat $MY_PIDFILE) " 
		if [ $? -eq 0 ]; then
			echo "$MY_NAME is running with PID "$($0 pid)
			exit 0
		else
			echo "$MY_NAME is not running (PIDFILE mismatch)"
			exit 1
		fi
		;;
	log|stdout)
		if [ -f $MY_LOGFILE ]; then
			tail -f $MY_LOGFILE
		else
			echo "No stdout output yet"
		fi
		;;
	
	err|stderr)
		if [ -f $MY_ERRFILE ]; then
			tail -f $MY_ERRFILE
		else
			echo "No stderr output yet"
		fi
		;;
	pid)
		if [ -f $MY_PIDFILE ]; then
			cat $MY_PIDFILE
		else
			echo "No pidfile found"
		fi
		;;
	run)
		daemonize
		payload
		;;
	help|?|--help|-h)
		echo "Usage: $0 [ start | stop | restart | (stat|status) | pause | resume | disable | enable | (log|stdout) | (err|stderr) ]"
		exit 0
		;;
	one)
		payload
		;;
	*)
		echo "Invalid argument"
		echo
		$0 help
		;;
esac
