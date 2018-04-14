#include "DongleInterface.h"

#include "hasp_hl.h"


/*#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


#ifdef _M_X64
#pragma comment(lib,"libhasp.lib")
#else
#pragma comment(lib,"libhasp_windows.lib")
#endif*/


CDongleInterface::CDongleInterface()
	: m_sVendorCode(
		"wVsGYIVDeUw/4sDL3Fr7mCrlIxxqeC85uaJ3GC7YipLzvOKCoZhIep95CJvm/ySbe5io48NHwrhIfvQl"
		"mxYteFzlS/YnOBN/dhKKbyJ73i2+w74smvlW95FryAzwDYH8G3bx7X3F1Eu3f1lNfU+UOsNCt19YwhS5"
		"Mdh+u6+aUx9I3meQvpVD+Po7QQty7FWKtvxkDemRf75yIBkYzcxmyXodcU6JgLXoUTc3OZjpR+ZY3fEo"
		"N12PFaAr0hwBPHcmPWavS+nSZc4ecVnsnJhOQ3kIkr2xhVP05flfWCD4cWsVdfRW4tZafxgb9qsbQXqi"
		"NPssjoiZwlnGcD1ZzDcTU7d+VwpWS9bXXtmEaUp7S5sYmNgeB4bEA/ssvDPW1HFHZFm3KlpfKYUrpYC8"
		"hmfGlcmuGyyiDBiN7cdeZDC7C0J8HYBagz7Kkv2x9m8wtrXAcs0eUVIT6Z+w+Ev2nQs06roskyJOLLu4"
		"S64kiwSihuPq/uG2iHMZ9mlORcVAu6eQ4KeddjVKKQxIXZpqYtktEtr6xq3Jw/EDCi06FRiuWbvhb0uX"
		"F5+88toVNCWUbsjZAtJqoFS00lME7INwtU0pRmchHm04EGFuuWKBp6B/5ebcdmqhAXqIY7NCitVgjNAW"
		"a+5RhAQM20v4JZF2Fl11Urjdj1KdoMh/S9U/OJZfgmFq1Y8ldJOnnbvu0o90C0VEOVvT6PPOM/W2SxEP"
		"MHI7/zviMz4uOw=="),
	  m_sErrorMessage(),
	  m_sDongleId()
{
}

CDongleInterface::~CDongleInterface()
{
}

bool CDongleInterface::ConnectDongle()
{
	bool connectSuccess(false);

	m_sErrorMessage = CString("");
	m_sDongleId = CString("Not Found");

	hasp_handle_t handle(0);
	const hasp_feature_t feature(HASP_PROGNUM_DEFAULT_FID | HASP_PROGNUM_OPT_NO_REMOTE);

	hasp_status_t status = hasp_login(feature, (hasp_vendor_code_t*)(LPCSTR)m_sVendorCode, &handle);

	/* check if operation was successful */
	switch (status)
	{
	case HASP_FEATURE_NOT_FOUND: // desired feature does not exist, or there is no licens in the dongle
	{
		m_sErrorMessage = "Feature does not exist, or licence expired";
	}
		break;
	case HASP_CONTAINER_NOT_FOUND: // vendor code was not found: dongle not present or wrong dongle
	{
		m_sErrorMessage = "Dongle not found or wrong dongle";
	}
		break;
	case HASP_OLD_DRIVER:
	{
		m_sErrorMessage = "Old driver found";
	}
		break;
	case HASP_NO_DRIVER:
	{
		m_sErrorMessage = "No driver found";
	}
		break;
	case HASP_INV_VCODE:
	{
		m_sErrorMessage = "Invalid vendor code";
	}
		break;
	case HASP_FEATURE_TYPE_NOT_IMPL:
	{
		m_sErrorMessage = "Requested feature not available";
	}
		break;
	case HASP_TMOF:
	{
		m_sErrorMessage = "Too many open handles";
	}
		break;
	case HASP_TS_DETECTED:
	{
		m_sErrorMessage = "Terminal server detected";
	}
		break;
	case HASP_INV_PROGNUM_OPT:
	{
		m_sErrorMessage = "Unknown program feature number requested";
	}
		break;
	case HASP_INSUF_MEM:
	{
		m_sErrorMessage = "Insufficient memory";
	}
		break;
	case HASP_STATUS_OK:
	{
		connectSuccess = true;

		char* info = 0;

		status = hasp_get_sessioninfo(handle, HASP_KEYINFO, &info);

		if (info)
		{	// extract dongle id
			char* start = strstr(info, "<haspid>");
			char* stop = strstr(info, "</haspid>");

			if (start && stop)
			{
				int idStart = (start + strlen("<haspid>")) - info;
				int idLen = (stop - info) - idStart;

				CStringA dongleId;

				char* data(dongleId.GetBufferSetLength(idLen));
				memcpy(data, info+idStart, idLen);

				m_sDongleId = dongleId;
			}

			hasp_free(info);
		}
	}
		break;
	default:
	{
		m_sErrorMessage = "Error reading dongle";
	}
		break;
	}

	hasp_logout(handle);

	return connectSuccess;
}

bool CDongleInterface::CheckFeature(int featNum)
{
	bool check(false);

	hasp_handle_t handle(0);
	const hasp_feature_t feature(featNum | HASP_PROGNUM_FEATURETYPE | HASP_PROGNUM_OPT_NO_REMOTE | HASP_PROGNUM_OPT_TS);

	hasp_status_t status = hasp_login(feature, (void*)(LPCSTR) m_sVendorCode, &handle);

	/* check if operation was successful */
	switch (status)
	{
	case HASP_FEATURE_NOT_FOUND: // desired feature does not exist, or there is no licens in the dongle
	{
		m_sErrorMessage = "Feature does not exist, or licence expired";
	}
		break;
	case HASP_CONTAINER_NOT_FOUND: // vendor code was not found: dongle not present or wrong dongle
	{
		m_sErrorMessage = "Dongle not found";
	}
		break;
	case HASP_OLD_DRIVER:
	{
		m_sErrorMessage = "Old driver found";
	}
		break;
	case HASP_NO_DRIVER:
	{
		m_sErrorMessage = "No driver found";
	}
		break;
	case HASP_INV_VCODE:
	{
		m_sErrorMessage = "Invalid vendor code";
	}
		break;
	case HASP_FEATURE_TYPE_NOT_IMPL:
	{
		m_sErrorMessage = "Requested feature not available";
	}
		break;
	case HASP_TMOF:
	{
		m_sErrorMessage = "Too many open handles";
	}
		break;
	case HASP_TS_DETECTED:
	{
		m_sErrorMessage = "Terminal server detected";
	}
		break;
	case HASP_INV_PROGNUM_OPT:
	{
		m_sErrorMessage = "Unknown program feature number requested";
	}
		break;
	case HASP_INSUF_MEM:
	{
		m_sErrorMessage = "Insufficient memory";
	}
		break;
	case HASP_STATUS_OK:
	{
		check = true;
	}
		break;
	default:
	{
		m_sErrorMessage = "Error reading dongle";
	}
		break;
	}

	hasp_logout(handle);

	return check;
}

CString CDongleInterface::GetErrorMessage()
{
	return m_sErrorMessage;
}

CString CDongleInterface::GetDongleId()
{
	return m_sDongleId;
}
