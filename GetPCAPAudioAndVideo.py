import os
import pyshark

# TODO: Find out what RTP payload type the VTC and VoIP guards will send in. They may be dynamically allocated.
PAYLOAD_TYPES = \
    {
        "audio" : [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19],
        "video" : [25, 26, 28, 31, 32, 34]
    }

class STREAM_TYPE:
    AUDIO = 0
    VIDEO = 1

def get_raw_audio(pcap_relative_name):
    rtp_payload_dict   = {}
    rtp_seq_dict       = {}
    ssrc_type          = {}
    protocol_filter    = "rtp"
    dirpath            = os.path.dirname(os.path.realpath(__file__))
    pcap_file          = f"{dirpath}/{pcap_relative_name}"

    print("Extracing rtp stream from: " + pcap_file)
    rtp_captures = pyshark.FileCapture(pcap_file, display_filter=protocol_filter)

    for rtp_capture in rtp_captures:
        rtp = rtp_capture[3]  # This element is a dictionary that contains the payload.

        if not hasattr(rtp, 'payload'):
            continue

        # Check for padding
        payload = rtp.payload.split(":")
        payload = payload if rtp.padding != 0 else payload[:-payload[-1]]

        # This implementation assumes the pcap has saved the audio in the correct sequence.
        if rtp.ssrc in rtp_payload_dict:
            rtp_payload_dict[str(rtp.ssrc)].append(payload)
            if int(rtp.seq) - 1 != rtp_seq_dict[str(rtp.ssrc)][-1]:
                print(f"Warning: RTP stream with source ID {rtp.ssrc} has packets missing.")
            rtp_seq_dict[str(rtp.ssrc)].append(int(rtp.seq))

        else:
            # Determine if stream is audio or video.
            p_type = int(rtp.p_type)
            if p_type not in PAYLOAD_TYPES["audio"] and p_type not in PAYLOAD_TYPES["video"]:
                raise NotImplementedError

            ssrc_type[str(rtp.ssrc)]        = STREAM_TYPE.AUDIO if p_type in PAYLOAD_TYPES["audio"] else STREAM_TYPE.VIDEO
            rtp_payload_dict[str(rtp.ssrc)] = [payload]
            rtp_seq_dict[str(rtp.ssrc)]     = [int(rtp.seq)]

    # Write raw bytes file for each ssrc in the pcap file.
    for rtp_ssrc, rtp_payload in rtp_payload_dict.items():
        outfile_name = "audio" if ssrc_type[rtp_ssrc] == STREAM_TYPE.AUDIO else "video"
        out_file     = f"{dirpath}/{outfile_name}_{rtp_ssrc}" + ".raw"

        with open(out_file,'wb') as raw_audio:
            for rtp_packet in rtp_payload:
                packet = " ".join(rtp_packet) # What is this bit, can I make it nicer?
                audio = bytearray.fromhex(packet)
                raw_audio.write(audio)

    print(f"\nFinished outputing in raw format.")

pcap_test = get_raw_audio("my.pcap")
