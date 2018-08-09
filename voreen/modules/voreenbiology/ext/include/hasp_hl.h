/*! \file hasp_hl.h HASP HL API declarations
 *
 */

/*!
 *  \mainpage HASP HL High Level API
 *  Copyright Aladdin Knowledge Systems Ltd.
 *
 *  $Id: hasp_hl.h,v 1.1.1.1 2008/01/11 16:08:40 Hockenbrink Exp $
 */


#ifndef __HASP_HL_H__
#define __HASP_HL_H__

#ifndef WITH_AKSTYPES
#if defined(_MSC_VER) || defined(__WATCOMC__) || defined(__BORLANDC__)
typedef unsigned __int64 hasp_u64_t;
typedef signed __int64 hasp_s64_t;
#else
typedef unsigned long long hasp_u64_t;
typedef signed long long hasp_s64_t;
#endif
#if defined(_MSC_VER)
typedef unsigned long hasp_u32_t;
typedef signed long hasp_s32_t;
#else
typedef unsigned int hasp_u32_t;
typedef signed int hasp_s32_t;
#endif
typedef unsigned short hasp_u16_t;
typedef signed short hasp_s16_t;
typedef unsigned char hasp_u8_t;
typedef signed char hasp_s8_t;
#endif

#if defined(WIN32) || defined(_MSC_VER) || defined(__BORLANDC__)
#define HASP_CALLCONV __stdcall
#else
#define HASP_CALLCONV
#endif

#ifdef __cplusplus
extern "C" {
#endif

/*! @defgroup hasp_general General declarations
 *
 * @{
 */

typedef hasp_u32_t hasp_status_t;        /*!< raw error code */
typedef hasp_u32_t hasp_size_t;          /*!< length */
typedef hasp_u32_t hasp_handle_t;        /*!< connection handle */
typedef hasp_u32_t hasp_feature_t;       /*!< feature id */
typedef hasp_u32_t hasp_fileid_t;        /*!< memory file id */
typedef hasp_u64_t hasp_time_t;          /*!< time, representing seconds since Jan-01-1970 0:00:00 GMT */

typedef void *hasp_vendor_code_t;       /*!< contains the vendor code */

#define HASP_UPDATEINFO  "<haspformat format=\"updateinfo\"/>"   /*!< hasp_get_sessioninfo() format to get update info (C2V) */
#define HASP_SESSIONINFO "<haspformat format=\"sessioninfo\"/>"  /*!< hasp_get_sessioninfo() format to get session info */
#define HASP_KEYINFO     "<haspformat format=\"keyinfo\"/>"      /*!< hasp_get_sessioninfo() format to get key/hardware info */

/*! @}
 */


/*! @defgroup hasp_feature_ids Feature ID defines
 *
 * See also \ref hasp_features
 *
 * @{
 */

/*! \brief "Featuretype" mask
 *
 *  AND-mask used to identify feature type
 */
#define HASP_FEATURETYPE_MASK 0xffff0000


/*! \brief "PROGRAM NUMBER FEATURE" type
 *
 *  After AND-ing with HASP_FEATURETYPE_MASK feature type contain this value.
 */
#define HASP_PROGNUM_FEATURETYPE 0xffff0000


/*! \brief program number mask
 *
 *  AND-mask used to extract program number from feature id if program number feature.
 */
#define HASP_PROGNUM_MASK 0x000000ff


/*! \brief prognum options mask
 *
 *  AND-mask used to identify prognum options:
 *     - HASP_PROGNUM_OPT_NO_LOCAL
 *     - HASP_PROGNUM_OPT_NO_REMOTE
 *     - HASP_PROGNUM_OPT_PROCESS
 *     - HASP_PROGNUM_OPT_CLASSIC
 *     - HASP_PROGNUM_OPT_TS
 *
 *  3 bits of the mask are reserved for future extensions and currently unused.
 *  Initialize them to zero.
 */
#define HASP_PROGNUM_OPT_MASK 0x0000ff00


/*! \brief "Prognum" option
 *
 *  Disable local license search
 */
#define HASP_PROGNUM_OPT_NO_LOCAL 0x00008000


/*! \brief "Prognum" option
 *
 *  Disable network license search
 */
#define HASP_PROGNUM_OPT_NO_REMOTE 0x00004000


/*! \brief "Prognum" option
 *
 *  Sets session count of network licenses to per-process
 */
#define HASP_PROGNUM_OPT_PROCESS 0x00002000


/*! \brief "Prognum" option
 *
 *  Enables the API to access "classic" (HASP4 or earlier) keys
 */
#define HASP_PROGNUM_OPT_CLASSIC 0x00001000


/*! \brief "Prognum" option
 *
 *  Presence of Terminal Services gets ignored
 */
#define HASP_PROGNUM_OPT_TS 0x00000800


/*! \brief HASP default feature id
 *
 *  Present in every hardware key.
 */
#define HASP_DEFAULT_FID 0


/*! \brief "Prognum" default feature id
 *
 *  Present in every hardware HASP key.
 */
#define HASP_PROGNUM_DEFAULT_FID (HASP_DEFAULT_FID | HASP_PROGNUM_FEATURETYPE)


/*! @}
 */

/*! \brief Minimal block size for hasp_encrypt() and hasp_decrypt() functions.
 */
#define HASP_MIN_BLOCK_SIZE 16

/*! \brief Minimal block size for legacy functions hasp_legacy_encrypt()
 *         and hasp_legacy_decrypt().
 */
#define HASP_MIN_BLOCK_SIZE_LEGACY 8

/*! @defgroup hasp_file_ids Memory file id defines
 *
 * @{
 */

/*! \brief HASP4 memory file
 *
 *  File id for HASP4 compatible memory contents w/o FAS
 */
#define HASP_FILEID_MAIN 0xfff0

/*! \brief HASP4 FAS memory file
 *
 *  (Dummy) file id for license data area of memory contents
 */
#define HASP_FILEID_LICENSE 0xfff2


/*! @}
 */

/*! @defgroup hasp_error_codes Error code defines
 *
 * @{
 */

enum hasp_error_codes
{
  HASP_STATUS_OK = 0,              /*!< no error occurred */
  HASP_MEM_RANGE = 1,              /*!< invalid memory address */
  HASP_INV_PROGNUM_OPT = 2,        /*!< unknown/invalid feature id option */
  HASP_INSUF_MEM = 3,              /*!< memory allocation failed */
  HASP_TMOF = 4,                   /*!< too many open features */
  HASP_ACCESS_DENIED = 5,          /*!< feature access denied */
  HASP_INCOMPAT_FEATURE = 6,       /*!< incompatible feature */
  HASP_CONTAINER_NOT_FOUND = 7,    /*!< license container not found */
  HASP_TOO_SHORT = 8,              /*!< en-/decryption length too short */
  HASP_INV_HND = 9,                /*!< invalid handle */
  HASP_INV_FILEID = 10,            /*!< invalid file id / memory descriptor */
  HASP_OLD_DRIVER = 11,            /*!< driver or support daemon version too old */
  HASP_NO_TIME = 12,               /*!< real time support not available */
  HASP_SYS_ERR = 13,               /*!< generic error from host system call */
  HASP_NO_DRIVER = 14,             /*!< hardware key driver not found */
  HASP_INV_FORMAT = 15,            /*!< unrecognized info format */
  HASP_REQ_NOT_SUPP = 16,          /*!< request not supported */
  HASP_INV_UPDATE_OBJ = 17,        /*!< invalid update object */
  HASP_KEYID_NOT_FOUND = 18,       /*!< key with requested id was not found */
  HASP_INV_UPDATE_DATA = 19,       /*!< update data consistency check failed */
  HASP_INV_UPDATE_NOTSUPP = 20,    /*!< update not supported by this key */
  HASP_INV_UPDATE_CNTR = 21,       /*!< update counter mismatch */
  HASP_INV_VCODE = 22,             /*!< invalid vendor code */
  HASP_ENC_NOT_SUPP = 23,          /*!< requested encryption algorithm not supported */
  HASP_INV_TIME = 24,              /*!< invalid date / time */
  HASP_NO_BATTERY_POWER = 25,      /*!< clock has no power */
  HASP_NO_ACK_SPACE = 26,          /*!< update requested acknowledgement, but no area to return it */
  HASP_TS_DETECTED = 27,           /*!< terminal services (remote terminal) detected */
  HASP_FEATURE_TYPE_NOT_IMPL = 28, /*!< feature type not implemented */
  HASP_UNKNOWN_ALG = 29,           /*!< unknown algorithm */
  HASP_INV_SIG = 30,               /*!< signature check failed */
  HASP_FEATURE_NOT_FOUND = 31,     /*!< feature not found */
  HASP_NO_LOG = 32,                /*!< trace log is not enabled */

  /* c++ use */
  HASP_INVALID_OBJECT = 500,
  HASP_INVALID_PARAMETER,
  HASP_ALREADY_LOGGED_IN,
  HASP_ALREADY_LOGGED_OUT,

  /* .net use */
  HASP_OPERATION_FAILED = 525,

  /* inside-api use */
  HASP_NO_EXTBLOCK = 600,          /*!< no classic memory extension block available */
  /* internal use */
  HASP_INV_PORT_TYPE = 650,        /*!< invalid port type */
  HASP_INV_PORT = 651,             /*!< invalid port value */
  /* catch-all */
  HASP_NOT_IMPL = 698,             /*!< capability isn't available */
  HASP_INT_ERR = 699               /*!< internal API error */
};

/*! @}
 */


/*! @defgroup hasp_basic The Basic API
 *
 * @{
 */

/* --------------------------------------------------------------------- */
/*! \brief Login into a feature.
 *
 *  This function establishes a context (logs into a feature).
 *
 *  \param feature_id       - Unique identifier of the feature\n
 *                            With "prognum" features (see \ref HASP_FEATURETYPE_MASK),
 *                            8 bits are reserved for legacy options (see \ref HASP_PROGNUM_OPT_MASK,
 *                            currently 5 bits are used):
 *                            - only local
 *                            - only remote
 *                            - login is counted per process ID
 *                            - disable terminal server check
 *                            - enable access to old (HASP3/HASP4) keys
 *  \param vendor_code      - pointer to the vendor code
 *  \param handle           - pointer to the resulting session handle
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_FEATURE_NOT_FOUND
 *                - the requested feature isn't available
 *              - HASP_CONTAINER_NOT_FOUND
 *                - no possible feature container found
 *              - HASP_FEATURE_TYPE_NOT_IMPL
 *                - the type of feature isn't implemented
 *              - HASP_INV_PROGNUM_OPT
 *                - unknown prognum option requested (\ref HASP_PROGNUM_OPT_MASK)
 *              - HASP_TMOF
 *                - too many open handles
 *              - HASP_INSUF_MEM
 *                - out of memory
 *              - HASP_INV_VCODE
 *                - invalid vendor code
 *              - HASP_NO_DRIVER
 *                - driver not installed
 *              - HASP_OLD_DRIVER
 *                - old driver installed
 *              - HASP_TS_DETECTED
 *                - program runs on a remote screen on Terminal Server
 *
 *  \sa  hasp_logout()
 *
 *  \remark
 *
 *  For local prognum features, concurrency is not handled and each login performs a decrement
 *  if it is a counting license.
 *
 *  Network prognum features just use the old HASPLM login logic with all drawbacks.
 *  There is only support for concurrent usage of \b one server (global server address).
 *
 */
hasp_status_t HASP_CALLCONV hasp_login(hasp_feature_t feature_id,
                                       hasp_vendor_code_t vendor_code,
                                       hasp_handle_t *handle);


/* --------------------------------------------------------------------- */
/*! \brief Logout.
 *
 *  Logs out from a session and frees all allocated resources for the session.
 *
 *  \param handle       - handle of session to log out from
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *
 *  \sa  hasp_login()
 */
hasp_status_t HASP_CALLCONV hasp_logout(hasp_handle_t handle);


/* --------------------------------------------------------------------- */
/*! \brief Encrypt a buffer.
 *
 *  This function encrypts a buffer.
 *
 *  \param handle      - session handle
 *  \param buffer      - pointer to the buffer to be encrypted
 *  \param length      - size in bytes of the buffer to be encrypted (16 bytes minimum)
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_TOO_SHORT
 *                - the length of the data to be encrypted is too short
 *              - HASP_ENC_NOT_SUPP
 *                - encryption type not supported by the hardware
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \sa hasp_decrypt()
 *
 *  \remark
 *  If the encryption fails (e.g. key removed in-between) the data pointed to by buffer is unmodified.
 */
hasp_status_t HASP_CALLCONV hasp_encrypt(hasp_handle_t handle, void *buffer, hasp_size_t length);


/* --------------------------------------------------------------------- */
/*! \brief Decrypt a buffer.
 *
 *  This function decrypts a buffer. This is the reverse operation of the
 *  hasp_encrypt() function. See \ref hasp_encrypt() for more information.
 *
 *  \param handle      - session handle
 *  \param buffer      - pointer to the buffer to be decrypted
 *  \param length      - size in bytes of the buffer to be decrypted (16 bytes minimum)
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_TOO_SHORT
 *                - the length of the data to be decrypted is too short
 *              - HASP_ENC_NOT_SUPP
 *                - encryption type not supported by the hardware
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \sa hasp_encrypt()
 *
 *  \remark
 *  If the decryption fails (e.g. key removed in-between) the data pointed to by buffer is unmodified.
 */
hasp_status_t HASP_CALLCONV hasp_decrypt(hasp_handle_t handle, void *buffer, hasp_size_t length);


/* --------------------------------------------------------------------- */
/*! \brief Read from key memory.
 *
 *  This function is used to read from the key memory.
 *
 *  \param      handle       - session handle
 *  \param      fileid       - id of the file to read (memory descriptor)
 *  \param      offset       - position in the file
 *  \param      length       - number of bytes to read
 *  \param      buffer       - result of the read operation
 *
 *  \return     status code.
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_INV_FILEID
 *                - unknown fileid
 *              - HASP_MEM_RANGE
 *                - attempt to read beyond eom
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \remark
 *  Valid fileids are \ref HASP_FILEID_LICENSE and \ref HASP_FILEID_MAIN.
 *
 *  \sa         hasp_write(), hasp_get_size()
 */
hasp_status_t HASP_CALLCONV hasp_read(hasp_handle_t handle, hasp_fileid_t fileid, hasp_size_t offset, hasp_size_t length, void *buffer);


/* --------------------------------------------------------------------- */
/*! \brief Write to key memory.
 *
 *  This function is used to write to the key memory. Depending on the provided
 *  session handle (either logged into the default feature or any other feature),
 *  write access to the FAS memory (\ref HASP_FILEID_LICENSE) is not permitted.
 *
 *  \param      handle       - session handle
 *  \param      fileid       - id of the file to write
 *  \param      offset       - position in the file
 *  \param      length       - number of bytes to write
 *  \param      buffer       - what to write
 *
 *  \return     status code.
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_INV_FILEID
 *                - unknown fileid
 *              - HASP_MEM_RANGE
 *                - attempt to read beyond eom
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \remark
 *  Valid fileids are \ref HASP_FILEID_LICENSE and \ref HASP_FILEID_MAIN.
 *
 *  \sa         hasp_read(), hasp_get_size()
 */
hasp_status_t HASP_CALLCONV hasp_write(hasp_handle_t handle, hasp_fileid_t fileid, hasp_size_t offset, hasp_size_t length, void *buffer);


/* --------------------------------------------------------------------- */
/*! \brief Get memory size.
 *
 *  This function is used to determine the memory size.
 *
 *  \param      handle       - session handle
 *  \param      fileid       - id of the file to query
 *  \param      size         - pointer to the resulting file size
 *
 *  \result     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_INV_FILEID
 *                - unknown fileid
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \remark
 *  Valid fileids are \ref HASP_FILEID_LICENSE and \ref HASP_FILEID_MAIN.
 *
 *  \sa         hasp_read(), hasp_write()
 */
hasp_status_t HASP_CALLCONV hasp_get_size(hasp_handle_t handle, hasp_fileid_t fileid, hasp_size_t *size);


/* --------------------------------------------------------------------- */
/*! \brief Read current time from a time key.
 *
 *  This function reads the current time from a time key.
 *  The time will be returned in seconds since Jan-01-1970 0:00:00 GMT.
 *
 *  \remark The general purpose of this function is not related to
 *          licensing, but to get reliable timestamps which are independent
 *          from the system clock.
 *  \remark This request is only supported on locally accessed keys. Trying to
 *          get the time from a remotely accessed key will return HASP_NO_TIME.
 *
 *  \param      handle       - session handle
 *  \param      time         - pointer to the actual time
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *              - HASP_NO_TIME
 *                - RTC not available or remote access
 *              - HASP_NO_BATTERY_POWER
 *                - RTC consistency check failed, battery probably dead
 *
 *  \sa hasp_datetime_to_hasptime(), hasp_hasptime_to_datetime()
 */
hasp_status_t HASP_CALLCONV hasp_get_rtc(hasp_handle_t handle, hasp_time_t *time);


/* --------------------------------------------------------------------- */
/*! @}
 */


/* --------------------------------------------------------------------- */
/*! @defgroup hasp_classic Legacy HASP functionality for backward compatibility
 *
 *  @{
 */

/*! \brief Legacy HASP4 compatible encryption function.
 *
 *  \param      handle       - session handle
 *  \param      buffer       - pointer to the buffer to be encrypted
 *  \param      length       - size in bytes of the buffer  (8 bytes minimum)
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_TOO_SHORT
 *                - the length of the data to be encrypted is too short
 *              - HASP_ENC_NOT_SUPP
 *                - encryption type not supported by the hardware
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \note The handle must have been obtained by calling \ref hasp_login() with
 *  a prognum feature id.
 *
 *  \sa hasp_legacy_decrypt(), hasp_encrypt(), hasp_decrypt()
 *
 *  \remark
 *  If the encryption fails (e.g. key removed in-between) the data pointed to by buffer is undefined.
 */
hasp_status_t HASP_CALLCONV hasp_legacy_encrypt(hasp_handle_t handle, void *buffer, hasp_size_t length);


/* --------------------------------------------------------------------- */
/*! \brief Legacy HASP4 compatible decryption function.
 *
 *  \param      handle       - session handle
 *  \param      buffer       - pointer to the buffer to be decrypted
 *  \param      length       - size in bytes of the buffer (8 bytes minimum)
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_TOO_SHORT
 *                - the length of the data to be decrypted is too short
 *              - HASP_ENC_NOT_SUPP
 *                - encryption type not supported by the hardware
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *
 *  \note The handle must have been obtained by calling \ref hasp_login() with
 *  a prognum feature id.
 *
 *  \sa hasp_legacy_encrypt(), hasp_decrypt(), hasp_encrypt()
 *
 *  \remark
 *  If the decryption fails (e.g. key removed in-between) the data pointed to by buffer is undefined.
 */
hasp_status_t HASP_CALLCONV hasp_legacy_decrypt(hasp_handle_t handle, void *buffer, hasp_size_t length);


/* --------------------------------------------------------------------- */
/*! \brief Write to HASP4 compatible real time clock
 *
 *  \param      handle       - session handle
 *  \param      new_time     - time value to be set
 *
 *  \remark This request is only supported on locally accessed keys. Trying to
 *          set the time on a remotely accessed key will return HASP_NO_TIME.
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *              - HASP_NO_TIME
 *                - RTC not available or access remote
 *
 *  \note The handle must have been obtained by calling \ref hasp_login() with
 *  a prognum feature id.
 *
 *  \sa hasp_get_rtc(), hasp_datetime_to_hasptime(), hasp_hasptime_to_datetime()
 */
hasp_status_t HASP_CALLCONV hasp_legacy_set_rtc(hasp_handle_t handle, hasp_time_t new_time);


/* --------------------------------------------------------------------- */
/*! \brief Set the LM idle time.
 *
 *  \param      handle       - session handle
 *  \param      idle_time    - the idle time in minutes
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *              - HASP_REQ_NOT_SUPP
 *                - attempt to set the idle time for a local license
 *
 *  \note The handle must have been obtained by calling \ref hasp_login() with
 *  a prognum feature id.
 */
hasp_status_t HASP_CALLCONV hasp_legacy_set_idletime(hasp_handle_t handle, hasp_u16_t idle_time);

/*! @}
 */


/*! @defgroup hasp_extended Extended HASP HL API
 *
 *  The extended API consists of functions which provide extended functionality. This
 *  advanced functionality is sometimes necessary, and addresses the "advanced" user.
 *
 *  @{
 */

/*! \brief Get information in a session context.
 *
 *  Memory for the information is allocated by this function and has to be freed by the
 *  \ref hasp_free() function.
 *
 *  \param      handle       - session handle
 *  \param      format       - XML definition of the output data structure
 *  \param      info         - pointer to the returned information (XML list)
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_CONTAINER_NOT_FOUND
 *                - the feature container isn't available anymore
 *              - HASP_INV_HND
 *                - invalid session handle
 *              - HASP_INV_FORMAT
 *                - unrecognized format
 *              - HASP_INSUF_MEM
 *                - out of memory
 *
 *  \sa hasp_free()
 *  \sa HASP_UPDATEINFO, HASP_SESSIONINFO, HASP_KEYINFO
 *  \sa \ref hasp_query_page
 */
hasp_status_t HASP_CALLCONV hasp_get_sessioninfo (hasp_handle_t handle, const char *format, char **info);


/*! \brief Free resources allocated by hasp_get_sessioninfo
 *
 *  This function must be called to free the resources allocated by hasp_get_sessioninfo()
 *  or to free the acknowledge data optionally returned from hasp_update().
 *
 *  \param      info   -  pointer to the resources to be freed
 *
 *  \sa hasp_get_sessioninfo(), hasp_update()
 */
void HASP_CALLCONV hasp_free (char *info);


/*! \brief Write an update.
 *
 *  This function writes update information. The update blob contains all necessary data
 *  to perform the update: Where to write (in which "container", e.g. dongle), the necessary
 *  access data (vendor code) and of course the update itself.
 *  If the update blob requested it, the function returns in an acknowledge blob,
 *  which is signed/encrypted by the updated instance and contains a proof that this update
 *  was successfully installed. Memory for the acknowledge blob is allocated by the API and has to be
 *  freed by the programmer (see \ref hasp_free()).
 *
 *  \param      update_data      - pointer to the complete update data.
 *  \param      ack_data         - pointer to a buffer to get the acknowledge data.
 *
 *  \return     status code
 *              - HASP_INV_UPDATE_DATA
 *                - required XML tags not found
 *                - contents in binary data missing or invalid
 *              - HASP_INV_UPDATE_OBJ
 *                - binary data doesn't contain an update blob
 *              - HASP_NO_ACK_SPACE
 *                - acknowledge data requested by the update, but ack_data input parameter is NULL
 *              - HASP_KEYID_NOT_FOUND
 *                - key to be updated not found
 *              - HASP_INV_UPDATE_NOTSUPP
 *                - update not supported by the key
 *              - HASP_INV_UPDATE_CNTR
 *                - update counter at the wrong position
 *              - HASP_INV_SIG
 *                - signature verification failed
 *
 *  \sa hasp_free()
 *
 *  \remark  Update via LM is not supported.
 */
hasp_status_t HASP_CALLCONV hasp_update(char *update_data, char **ack_data);

/*! @}
 */


/* --------------------------------------------------------------------- */
/*! @defgroup hasp_util Utility functions
 *
 *  @{
 */

/* --------------------------------------------------------------------- */
/*! \brief Convert broken up time into a time type
 *
 *  \param      day          - input day
 *  \param      month        - input month
 *  \param      year         - input year
 *  \param      hour         - input hour
 *  \param      minute       - input minute
 *  \param      second       - input second
 *  \param      time         - pointer to put result
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_TIME
 *                - time outside of the supported range
 *
 *  \remark Times are in UTC.
 *
 *  \sa hasp_hasptime_to_datetime()
 */
hasp_status_t HASP_CALLCONV hasp_datetime_to_hasptime(unsigned int day,
                                                      unsigned int month,
                                                      unsigned int year,
                                                      unsigned int hour,
                                                      unsigned int minute,
                                                      unsigned int second,
                                                      hasp_time_t *time);


/* --------------------------------------------------------------------- */
/*! \brief Convert time type into broken up time
 *
 *  \param      time         - pointer to put result
 *  \param      day          - pointer to day
 *  \param      month        - pointer to month
 *  \param      year         - pointer to year
 *  \param      hour         - pointer to hour
 *  \param      minute       - pointer to minute
 *  \param      second       - pointer to second
 *
 *  \return     status code
 *              - HASP_STATUS_OK
 *                - the request completed successfully
 *              - HASP_INV_TIME
 *                - time outside of the supported range
 *
 *  \remark Times are in UTC.
 *
 *  \sa hasp_datetime_to_hasptime()
 */
hasp_status_t HASP_CALLCONV hasp_hasptime_to_datetime(hasp_time_t time,
                                                      unsigned int *day,
                                                      unsigned int *month,
                                                      unsigned int *year,
                                                      unsigned int *hour,
                                                      unsigned int *minute,
                                                      unsigned int *second);


/*! @}
 */


/*! \page hasp_features Feature ID convention

\section FeatIntr Feature ID  introduction

Feature ids are 32bits wide. If the upper 16 bit contain the value indicated by \ref HASP_PROGNUM_FEATURETYPE,
the feature defines a prognum feature.\n
For prognum features there are some options encoded in the feature id. These include
  - \ref HASP_PROGNUM_OPT_NO_LOCAL\n
    Don't search for a license locally. "Remote-only"\n\n
  - \ref HASP_PROGNUM_OPT_NO_REMOTE\n
    Don't search for a license on the network. "Local-only"\n\n
  - \ref HASP_PROGNUM_OPT_PROCESS\n
    In case, a license is found in the network, count license usage per process
    instead of per workstation.\n\n
  - \ref HASP_PROGNUM_OPT_TS\n
    Don't detect whether the program is running on a remote screen on a terminal server.\n\n
  - \ref HASP_PROGNUM_OPT_CLASSIC\n
    The API by default only searches for HASPHL keys. When this option is set, it also
    searches for HASP3/HASP4 keys.\n\n
*/


/*! \page hasp_query_page hasp_get_sessioninfo format and info

\section QueryPage hasp_get_sessioninfo format and info

Calling hasp_get_sessioninfo() with \ref HASP_UPDATEINFO format will return something like this:

\verbatim
<?xml version="1.0" encoding="UTF-8"?>
<hasp_info>
        <c2v>
                YYIBlIADY3R2oQaABEAwulCiCYABaoEBBIIBAKOBxoABAIGBwD2sfFj8UKuDvNWH9
                LhfRKDzUbLCAi6E9mN8ea7EclwOl9VeLMDuLvfsEvkor2igmwxg/wWs6HCuypEFi6
                V/FkI4EUmQNmcKSIY302s9CzHP7aCrG7QKvzArVq25Nc7UxIQJ4kZJm1oWiw3zZJq
                UY+G0EleETkPZ8n2uDfMauBpdWhW0R35rHlRM4wiYCZzaelpRtDX36HDh1caqfpaL
                mUnwWXRz0+tLs+Dvd+kLmvcQ6jWJJb4r2rxywG2IW1WTjIWBsI+h0/UgaIhG1J+9R
                EQ1SrMx3YQ2bpdlK3FluZVDayW9okv7idxKJS4zGG+4UOccpKT4aWJi9cR0vdm4s/
                J6fUNbhK522x/gdvR51a6ll46GpVn2HjD0ZpAgCeu6xAIwHJ7Kc6tjeRfxYX9YksE
                aB9JoV/uaPTHnbu2AgQmd0r09p0zmXgD4Kuk8EtTs1GoBbY7WF3qHJsj1Iz1ZeAdA
                rdNOYKsOgA/q1tuLLR7O0dag
        </c2v>
</hasp_info>
\endverbatim


Calling hasp_get_sessioninfo() with \ref HASP_SESSIONINFO format will return something like this:

\verbatim
<hasp_info>
        <feature>
                <featureid>4294905856</featureid>
                <maxlogins>5</maxlogins>
                <currentlogins>1</currentlogins>
                <activations>unlimited</activations>
        </feature>
</hasp_info>
\endverbatim

In case of a expiring license on a time enabled key (prognum <= 8),
instead of the remaining activations the expiration date will be
returned:

\verbatim
<hasp_info>
        <feature>
                <featureid>4294905857</featureid>
                <expirationdate>1052919239</expirationdate>
        </feature>
</hasp_info>
\endverbatim

For locally accessed keys there is no \p maxlogins and \p currentlogins field.


Calling hasp_get_sessioninfo() with \ref HASP_KEYINFO format will return something like
this for a locally accessed key:

\verbatim
<?xml version="1.0" encoding="UTF-8"?>
<hasp_info>
        <keyspec>
                <keycaps>
                        <hasp4/>
                        <aes/>
                        <newintf/>
                </keycaps>
                <hasp>
                        <haspid>12345</haspid>
                        <nethasptype>0</nethasptype>
                        <memoryinfo>
                                <name>"Main"</name>
                                <fileid>65520</fileid>
                                <size>48</size>
                        </memoryinfo>
                        <memoryinfo>
                                <name>"FAS"</name>
                                <fileid>65522</fileid>
                                <size>80</size>
                        </memoryinfo>
                </hasp>
                <port>
                        <type>"USB"</type>
                        <address>1</address>
                </port>
        </keyspec>
</hasp_info>
\endverbatim


Calling hasp_get_sessioninfo() with \ref HASP_KEYINFO format will return something like
this for a remotely accessed key:

\verbatim
<?xml version="1.0" encoding="UTF-8"?>
<hasp_info>
        <keyspec>
                <keycaps>
                        <hasp4/>
                </keycaps>
                <hasp>
                        <haspid>782062012</haspid>
                        <nethasptype>5</nethasptype>
                        <memoryinfo>
                                <name>"Main"</name>
                                <fileid>65520</fileid>
                                <size>432</size>
                        </memoryinfo>
                        <memoryinfo>
                                <name>"FAS"</name>
                                <fileid>65522</fileid>
                                <size>448</size>
                        </memoryinfo>
                </hasp>
                <serveraddress>
                        <protocol>"IP"</protocol>
                        <address>"10.20.3.10"</address>
                </serveraddress>
        </keyspec>
</hasp_info>
\endverbatim

\p keycaps flags:
         - hasp4
           - support HASP4 compatible encryption
         - aes
           - support AES encryption
         - rtc
           - key has real time clock chip
         - newintf
           - supports new access interface

*/

#ifdef __cplusplus
}
#endif

#endif /* #ifndef __HASP_HL_H__ */
