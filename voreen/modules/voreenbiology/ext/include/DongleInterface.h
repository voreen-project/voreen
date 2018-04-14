#ifndef VRN_DONGLEINTERFACE_H
#define VRN_DONGLEINTERFACE_H

#include "voreen/core/voreencoreapi.h"

#include "atlstr.h"

/**
 * The class CDongleInterface communicates with the dongle. 
 * It sets up a connection to the dongle and checks the feature present in the dongle
 */
class VRN_CORE_API CDongleInterface
{
public:
	CDongleInterface();
	~CDongleInterface();

	/**
	 * Connects to the dongle and checks for the venor code
	 */
	bool ConnectDongle();

	/**
	 * Check, if the given feature is active in the dongle
	 */
	bool CheckFeature(int featNum);

	/**
	 * Return last error message
	 */
	CString GetErrorMessage();

	/**
	 * Return dongle id
	 */
	CString GetDongleId();

private:
	CStringA m_sVendorCode;

	CString m_sErrorMessage;
	CString m_sDongleId;
};

#endif
