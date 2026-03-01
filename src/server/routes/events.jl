# © 2026 Joshua Benjamin Jewell. All rights reserved.
# Licensed under the GNU Affero General Public License version 3 (AGPLv3).

# Routes: /api/v1/events  (global SSE stream)

using HTTP

@get "/api/v1/events" function(req)
    ch = subscribe_events()
    headers = [
        "Content-Type"  => "text/event-stream",
        "Cache-Control" => "no-cache",
        "Connection"    => "keep-alive",
        "X-Accel-Buffering" => "no",
    ]
    # Oxygen streams the response body from a Channel or IO.
    # We open a pipe and write SSE frames as events arrive.
    body = IOBuffer()
    try
        # Initial heartbeat so the client knows the connection is live
        write(body, ": connected\n\n")
        while isopen(ch)
            if isready(ch)
                write(body, take!(ch))
            else
                # Keepalive comment every 15s prevents proxy timeouts
                write(body, ": keepalive\n\n")
                sleep(15)
            end
        end
    finally
        unsubscribe_events(ch)
    end
    HTTP.Response(200, headers; body=take!(body))
end
